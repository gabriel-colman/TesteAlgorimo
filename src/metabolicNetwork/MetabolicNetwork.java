/* 
PITUFO - A software tool to find all minimal precursor sets for a given set of targets in metabolic networks.

Copyright (C) 2011 Ludovic Cottret (l.cottret@gmail.com <mailto:l.cottret@gmail.com>), Paulo Vieira Milreu  (paulovieira@milreu.com.br <mailto:paulovieira@milreu.com.br>)      
This file is part of PITUFO.

PITUFO is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PITUFO is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PITUFO.  If not, see <http://www.gnu.org/licenses/>.
 */
package metabolicNetwork;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jdom2.Document;
import org.jdom2.Element;
import org.jdom2.JDOMException;
import org.jdom2.input.SAXBuilder;
import org.sbml.jsbml.Model;
import org.sbml.jsbml.SBMLDocument;
import org.sbml.jsbml.SBMLReader;


//import pitufo.ParallelFeasibilityTest;
import application.InputParameters;

public class MetabolicNetwork {

	// Compounds and reactions definition    
	HashMap<String, Compound> compounds = new HashMap<String, Compound>();
	HashMap<String,Reaction>  reactions = new HashMap<String, Reaction>();

	public MetabolicNetwork hardCopy() {

		MetabolicNetwork copy = new MetabolicNetwork();
		// copy the compounds
		for(Compound c: compounds.values()) {
			if(c.hasHighDegree() == true){
				System.out.println(c.getId() + " had a high degree ");
			}
			Compound copiedCompound = copy.addCompound(c.getId(), c.getName(), c.compartment);
			copiedCompound.copyPropertiesFrom(c);
		}

		// copy the reactions
		int nbReactions = 0;
		for(Reaction r: reactions.values()) {
			if( r.isReversible() && r.getId().endsWith("_REV") )
				continue;

			Reaction copiedReaction = copy.addNewReaction(r.getId(), r.getName(), r.isReversible());
			copiedReaction.setIdxBitSet(nbReactions++);
			copiedReaction.addCofactors(r.getCofactors());
			Reaction r2 = null;
			if( r.isReversible() ) {
				r2 = copy.addNewReaction(r.getId()+"_REV", r.getName()+"_REV", r.isReversible());
				r2.setIdxBitSet(nbReactions++);
				copiedReaction.setReverse(r2);
				r2.addCofactors(r.getCofactors());
			}

			// add the substrates of the reaction
			for(Compound substrate: r.getSubstrates().values()) {
				// Finds the compound
				String cId = substrate.getId();
				Compound c = copy.getCompounds().get( cId );
				if( c != null)
				{
					Double stoich = 1.0;
					if(r.getSubstratesStoich().containsKey(cId)){
						stoich = r.getSubstratesStoich().get(cId);
					}

					copiedReaction.addSubstrate(c, stoich);
					if( r.isReversible() ){
						r2.addProduct(c, stoich);
					}
				}
			}

			// add the products of the reaction
			for(Compound product: r.getProduces().values()) {
				// Finds the compound
				String cId = product.getId();
				Compound c = copy.getCompounds().get( cId );
				if( c != null)
				{
					Double stoich = 1.0;
					if(r.getProductsStoich().containsKey(cId)){
						stoich = r.getProductsStoich().get(cId);
					}
					copiedReaction.addProduct(c, stoich);
					if( r.isReversible() ){
						r2.addSubstrate(c, stoich);
					}
				}
			}        
		}

		//update the consumedBy and producedBy lists for each compound
		for(Reaction r : copy.getReactions().values()){
			for(Compound c : r.getSubstrates().values()){
				c.addSubstrateOf(r);
			}
			for(Compound c : r.getProduces().values()){
				c.addProducedBy(r);
			}
		}

		return copy;
	}


	public MetabolicNetwork backupHardCopy() {

		MetabolicNetwork copy = new MetabolicNetwork();
		// copy the compounds
		for(Compound c: compounds.values()) {
			//Compound copiedCompound = new Compound(c.getId(), c.getName(), c.compartment);
			Compound copiedCompound = copy.addCompound(c.getId(), c.getName(), c.compartment);
			copiedCompound.copyPropertiesFrom(c);

		}

		// copy the reactions
		for(Reaction r: reactions.values()) {
			if( r.isReversible() && r.getId().endsWith("_REV") )
				continue;

			Reaction copiedReaction = copy.addNewReaction(r.getId(), r.getName(), r.isReversible());
			Reaction r2 = null;
			if( r.isReversible() ) {
				r2 = copy.addNewReaction(r.getId()+"_REV", r.getName()+"_REV", r.isReversible());
				copiedReaction.setReverse(r2);
			}

			// add the substrates of the reaction
			for(Compound substrate: r.getSubstrates().values()) {
				// Finds the compound
				String cId = substrate.getId();
				Compound c = copy.getCompounds().get( cId );
				if( c != null)
				{
					Double stoich = 1.0;
					if(r.getSubstratesStoich().containsKey(cId)){
						stoich = r.getSubstratesStoich().get(cId);
					}
					copiedReaction.addSubstrate(c, stoich);
					if( r.isReversible() )
						r2.addProduct(c, stoich);
				}
			}

			// add the products of the reaction
			for(Compound product: r.getProduces().values()) {
				// Finds the compound
				String cId = product.getId();
				Compound c = copy.getCompounds().get( cId );
				if( c != null)
				{
					Double stoich = 1.0;
					if(r.getProductsStoich().containsKey(cId)){
						stoich = r.getProductsStoich().get(cId);
					}
					copiedReaction.addProduct(c, stoich);
					if( r.isReversible() )
						r2.addSubstrate(c, stoich);
				}
			}        
		}

		return copy;
	}

	public void print(){
		for (Iterator<Compound> iter = compounds.values().iterator(); iter.hasNext();) 
		{
			Compound c = (Compound) iter.next();
			System.out.println(c.id);

			System.out.print("Substrate->");
			for (Iterator<Reaction> iterator = c.getSubstrateOf().iterator(); iterator.hasNext();) 
			{
				Reaction r = (Reaction) iterator.next();
				System.out.print(r.id + " ");
			}

			System.out.print("\nProduced By->");
			for (Iterator<Reaction> iterator = c.getProducedBy().iterator(); iterator.hasNext();) 
			{
				Reaction r = (Reaction) iterator.next();
				System.out.print(r.id + " ");
			}
			System.out.println("\nTarget: "+c.isTarget());
			System.out.println("Precursor: "+c.isPrecursor());
			System.out.println("Bootstrap: "+c.isBootstrap());

			System.out.print("\n");
		}

		for (Iterator<Reaction> iter = reactions.values().iterator(); iter.hasNext();) 
		{
			Reaction r = (Reaction) iter.next();
			System.out.println(r.id);
		}
	}

	public void readFromSbmlFormat(String inputFile)
	{
		try
		{
			SBMLDocument xmlNetwork = (new SBMLReader()).readSBML(inputFile);
			Model network = xmlNetwork.getModel();

			addCompounds(network);
			addReactions(network);
		}
		catch(Exception e){
			System.out.println(e.getMessage());
		}

	}
/*
	private void filterPairedCofactors(String inputFile) {
		// TODO Auto-generated method stub
		BufferedReader br = null;

		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(inputFile));

			boolean foundCofactors = false;
			String reactionID = "";
			List<String> coFactors = new ArrayList<String>();
			while ((sCurrentLine = br.readLine()) != null) {
				if(sCurrentLine.contains("<reaction id")){
					reactionID = sCurrentLine.replaceAll("(.*<reaction id=\")(.+?)(\".*)", "$2");
					System.out.println("Reac: " + reactionID);
				}
				else if(sCurrentLine.contains("<cofactors>")){
					foundCofactors = true;
				}
				else if(sCurrentLine.contains("</cofactors>")){
					// find reaction object
					Reaction r = this.reactions.get(reactionID);
					for(String reac : this.reactions.keySet()){
						System.out.println("Search for " + reac);
					}
					if(r == null){
						System.exit(0);
					}

					// check if C atoms are still balanced when we exclude co-factors
					int nbCAtomsInSubstrates = 0;
					int nbCAtomsInProducts = 0;
					List<Compound> coFactorCompoundsSubs = new ArrayList<Compound>();
					List<Compound> coFactorCompoundsProd = new ArrayList<Compound>();
					for(Compound c : r.getSubstrates().values()){
						if(! coFactors.contains(c.getId())){ // not a co-factor
							nbCAtomsInSubstrates += c.getNumberOfAtom("C") * r.getSubstrateStochiometricValue(c);
						}
						else{
							coFactorCompoundsSubs.add(c);
						}
					}
					for(Compound c : r.getProduces().values()){
						if(! coFactors.contains(c.getId())){ // not a co-factor
							nbCAtomsInProducts += c.getNumberOfAtom("C") * r.getProductStochiometricValue(c);
						}
						else{
							coFactorCompoundsProd.add(c);
						}
					}

					// still balanced
					if(nbCAtomsInSubstrates == nbCAtomsInProducts){
						System.out.println("Reaction " + r + " is balanced");
						for(Compound c : coFactorCompoundsSubs){
							r.removeSubstrate(c);
							System.out.println("Remove " + c.getId());
						}
						for(Compound c : coFactorCompoundsProd){
							r.removeProduct(c);
							System.out.println("Remove " + c.getId());
						}

						// are there still substrates and products in the reaction ?
						if(r.getSubstrates().size() == 0 || r.getProduces().size() == 0){
							this.reactions.remove(r.getId());
							System.out.println("Removed reaction because there are no subtrates/products left");
						}
					}
					else{
						System.out.println("Reaction " + r + " is not balanced");
					}


					foundCofactors = false;
					coFactors.clear();
				}
				else if(foundCofactors == true && sCurrentLine.contains("<speciesReference species")){
					String coFactor = sCurrentLine.replaceAll("(.*<speciesReference species=\")(.+?)(\".*)", "$2");
					System.out.println("CoFactor: " + coFactor);
					coFactors.add(coFactor);
				}
			}

		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
	}
*/

	/*
	private void getSmiles(String inputFile) {
		// TODO Auto-generated method stub

		BufferedReader br = null;

		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(inputFile));

			while ((sCurrentLine = br.readLine()) != null) {
				if(sCurrentLine.contains("<species id") && sCurrentLine.contains("formula")){
					//get smiles formula
					String formula = sCurrentLine.replaceAll("(.*formula=\")(.*?)(\".*)", "$2");
					if(! formula.equals("")){ // -> there is a formula
						// get compound ID
						String id = sCurrentLine.replaceAll("(.*<species id=\")(.*?)(\".*)", "$2");
						System.out.println("ID " + id + " -> " + formula);

						// get atoms and number of atoms 
						String[] atoms = formula.split("\\d+"); // split on numbers
						String[] nbAtoms = formula.split("\\D+"); // split on letters

						// build hash map AtomName => #atoms
						HashMap<String, Integer> atomMap = new HashMap<String, Integer>();
						for(int i = 0; i < atoms.length; i++){
							atomMap.put(atoms[i], Integer.parseInt(nbAtoms[i+1])); // nbAtoms[i+1], because nbAtoms[0] is empty -> string before formula
						}
						System.out.println(atomMap);

						// find the compound object that corresponds to the id
						boolean foundID = false;
						for(String cID: this.compounds.keySet()){
							if(cID.equals(id)){
								foundID = true;
								// add atoms to the compound object
								this.compounds.get(cID).setAtoms(atomMap);
								break;
							}
						}
						if(foundID == false){
							System.out.println("Could not find compound " + id + " to add the formula");
							System.exit(0);
						}
					}
				}
			}

		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
	}
*/

	public Compound addCompound(String id, String name, String compartment)
	{
		Compound c = new Compound(id, name, compartment);
		compounds.put(id, c);
		return c;
	}

	public void addCompounds(Model model)
	{
		int nc = model.getNumSpecies();
		for(int i = 0; i < nc; i++) {
			addCompound(model.getSpecies(i).getId(), model.getSpecies(i).getName(), model.getSpecies(i).getCompartment());
		}
	}

	public Reaction addNewReaction(String id, String name, boolean reversible)
	{
		Reaction r = new Reaction(id, name, reversible);
		if(reactions.containsKey(id)){
			System.out.println("Already in hash " + id);
		}
		reactions.put(id, r);
		return r;	
	}

	public void removeReactionDependencies(Reaction r)
	{
		for(Compound c: r.produces.values())
			c.removeProducedBy(r);
		//c.getProducedBy().remove(r);
		for(Compound c: r.substrates.values()){
			c.removeSubstrateOf(r);
		}
		//c.getSubstrateOf().remove(r);
	}

	public void removeReactionAndItsReverse(Reaction r){
		removeReactionDependencies(r);
		if( r.isReversible() && r.getReverseReaction() != null ) {
			removeReactionDependencies(r.getReverseReaction());
		}

		/*if(r.isReversible()){
			reactions.remove(r.getReverseReaction());
		}
		reactions.remove(r);
		 */

		Iterator<Reaction> iter = reactions.values().iterator();
		while(iter.hasNext()){
			Reaction nextReac = iter.next();
			if(nextReac.getId().equals(r.getId()) || (r.isReversible() && nextReac.getId().equals(r.getId() + "_REV"))){
				iter.remove();
				System.out.println("Remove " + nextReac + " from network");
			}
		}
	}

	public void removeReaction(Reaction r)
	{
		removeReactionDependencies(r);
		if( r.isReversible() && r.getReverseReaction() != null ) {
			r.getReverseReaction().setReversible(false);
		}
		//reactions.values().remove(r);
		Iterator<Reaction> iter = reactions.values().iterator();
		while(iter.hasNext()){
			if(iter.next().getId().equals(r.getId())){
				iter.remove();
				break;
			}
		}
	}



	public void recursivelyRemoveReaction(Reaction r)
	{
		List<Compound> products = new ArrayList<Compound>(r.getProduces().values());
		removeReactionDependencies(r);
		if( r.isReversible() && r.getReverseReaction() != null ) {
			r.getReverseReaction().setReversible(false);
		}
		Iterator<Reaction> iter = reactions.values().iterator();
		while(iter.hasNext()){
			if(iter.next().getId().equals(r.getId())){
				iter.remove();
				break;
			}
		}
		for (Compound c : products){
			Compound realC = this.compounds.get(c.getId());
			if (!realC.isPrecursor() && realC.getProducedBy().size() == 0){
				List<Reaction> reactionToRemove = new ArrayList<Reaction>();
				for (Reaction nr : c.getSubstrateOf()){
					reactionToRemove.add(nr);
				}
				for(Reaction rr : reactionToRemove){
					recursivelyRemoveReaction(rr);
				}
			}
		}
	}

	public void addReactions(Model network)
	{
		int nr = network.getNumReactions();
		for(int i = 0; i < nr; i++) {
			Reaction r1 = addNewReaction(network.getReaction(i).getId(), network.getReaction(i).getName(), network.getReaction(i).isReversible());
			Reaction r2 = null;
			if( r1.isReversible() ) {
				r2 = addNewReaction(r1.getId()+"_REV", r1.getName()+"_REV", r1.isReversible());
				r1.setReverse(r2);
			}

			// add the substrates of the reaction
			int ns = network.getReaction(i).getNumReactants();
			for(int j = 0; j < ns; j++) {
				// Finds the compound
				String cId = network.getReaction(i).getReactant(j).getSpecies();
				double stoich = network.getReaction(i).getReactant(j).getStoichiometry();
				Compound c = compounds.get( cId );
				if( c != null)
				{
					//r1.addSubstrate(c);
					r1.addSubstrate(c, stoich);
					if( r1.isReversible() )
						r2.addProduct(c, stoich);
				}

			}

			// add the products of the reaction
			int np = network.getReaction(i).getNumProducts();
			for(int j = 0; j < np; j++) {
				// Finds the compound
				String cId = network.getReaction(i).getProduct(j).getSpecies();
				double stoich = network.getReaction(i).getProduct(j).getStoichiometry();
				Compound c = compounds.get(cId);
				if( c != null)
				{
					r1.addProduct(c, stoich);
					if( r1.isReversible() )
						r2.addSubstrate(c, stoich);
				}
			}        
		}
	}

	/*
	 * Returns true if the network can produce the compound c and false otherwise.
	 * 
	 */
	public boolean canProduce(Compound c)
	{
		for(Reaction r: reactions.values())
		{
			if( r.canSynthetize(c)  )
				return true;
		}
		return false;
	}


	/*
	 * Removes the compound c. If "removeReaction" is true, the reactions involving c are also removed. 
	 * 
	 */
	public void removeCompound(Compound c, boolean removeReaction)
	{
		if( removeReaction )
		{
			List<Reaction> reactionsInvolvingC = new ArrayList<Reaction>();
			for(Reaction r: reactions.values())
			{
				if( r.substrates.containsKey(c.getId()) || r.produces.containsKey(c.getId()) )
					reactionsInvolvingC.add(r);
			}
			for(int i = reactionsInvolvingC.size()-1; i >= 0; i--)
				removeReaction(reactionsInvolvingC.get(i));

		}
		compounds.values().remove(c);
	}

	/**
	 * Removes compounds from the network
	 * @param cpdsToRemove
	 */
	public void removeCompounds(List<Compound> cpdsToRemove) {

		HashMap<String, Compound> newListOfCompounds = new HashMap<String, Compound>();
		for(Compound cpd : compounds.values()) 
		{
			if(!cpdsToRemove.contains(cpd)) 
				newListOfCompounds.put(cpd.getId(), cpd);
		}
		compounds = newListOfCompounds;

		// update reactions removing the deleted compounds from their product/substrate lists.
		for(Reaction reaction : reactions.values()) {
			HashMap<String, Compound> newSubstrates = new HashMap<String, Compound>();
			for(Compound cpd : reaction.getSubstrates().values()) 
			{
				if(compounds.containsKey(cpd.getId())) {
					newSubstrates.put(cpd.getId(), cpd);
				}
			}
			reaction.substrates = newSubstrates;

			HashMap<String, Compound> newProducts = new HashMap<String, Compound>();
			for(Compound cpd : reaction.getProduces().values()) 
			{
				if(compounds.containsKey(cpd.getId())) {
					newProducts.put(cpd.getId(), cpd);
				}
			}
			reaction.produces = newProducts;
		}
		//cleanEmptyReactions();
	}

	// Looks for compounds "c" marked as "non topological precursors" and transforms them into topological precursors
	// by adding a fake reaction "c_precursor -> c" and no reaction that produces "c_precursor".
	public void extendNonTopologicalPrecursors()
	{
		if (!InputParameters.getConsiderOnlyUserDefinedPrecursors()){
			// Identify all the non topological precursors in the network
			List<Compound> nonTopologicalPrecursors = new ArrayList<Compound>();
			for(Compound c : getCompounds().values() )
			{
				if( c.isPrecursor() && !c.isTopologicalPrecursor() )
					nonTopologicalPrecursors.add(c);
			}

			// And adds to them the special reactions
			for(Compound c : nonTopologicalPrecursors )
			{
				// Create an special reaction that takes as substrate a special compound c_precursor  
				// and that produces the compound c.
				Compound precursor =addCompound("PRECURSOR_"+c.getId(), "PRECURSOR_"+c.getId(), "PRECURSOR_"+c.getId());
				precursor.setTopologicalPrecursor(true);
				Reaction special = addNewReaction("PRECURSOR_"+c.getId(), "PRECURSOR_"+c.getId(), false);
				special.addProduct(c);
				special.addSubstrate(precursor);
			}
		}
	}	

	public void extendTopologicalPrecursors(boolean precursorIfProducedOnlyByReversible) 
	{
		if (!InputParameters.getConsiderOnlyUserDefinedPrecursors()){
			for(Compound cpd : this.compounds.values())
			{
				List<Reaction> reactionsThatProduce = cpd.getReactionsThatProduce(false);
				if( reactionsThatProduce.size() == 0 )
				{
					cpd.setTopologicalPrecursor(true);
				}
				// not sure if it is a precursor when there is only one reaction that consumes
				else if( precursorIfProducedOnlyByReversible && reactionsThatProduce.size() == 1 && cpd.getReactionsThatConsume(false).size() == 1 && reactionsThatProduce.get(0).isReversible() )
					//else if( precursorIfProducedOnlyByReversible && reactionsThatProduce.size() == 1 && reactionsThatProduce.get(0).isReversible() )
				{
					cpd.setTopologicalPrecursor(true);
				}
			}
		}
	}

	public void removeTopologicalSourcesThatAreNotPrecursors() 
	{
		boolean changed = true;
		List<Compound> cpdsToRemove = new ArrayList<Compound>();
		while( changed ) {
			changed = false;
			cpdsToRemove.clear();
			for(Compound cpd : this.compounds.values()) 
			{
				List<Reaction> reactionsThatProduce = cpd.getReactionsThatProduce(true);
				if( reactionsThatProduce.size() == 0 && !cpd.isPrecursor() )
				{
					cpdsToRemove.add(cpd);
				}
			}

			if(cpdsToRemove.size() > 0) {
				changed = true;
				for(Compound c: cpdsToRemove) {
					removeCompound(c, true);
				}
			}
		}
	}

	/*
	 * This procedure eliminates reactions which has product or substrates side empty
	 */
	public void cleanEmptyReactions() {
		List<Reaction> reactionsToRemove = new ArrayList<Reaction>();

		// for each reaction checks whether their products or substrates are empty
		for(Reaction r: getReactions().values()) {
				if( r.getProduces().isEmpty()) {
					reactionsToRemove.add(r);
				}
				if (r.getSubstrates().isEmpty()){
					reactionsToRemove.add(r);
				}
		}

		for(Reaction rToDel: reactionsToRemove) {
			System.out.println("Without subs/prod " + rToDel);
			removeReaction(rToDel);
			//removeReactionAndItsReverse(rToDel);
		}		
	}	

	/*
	 * This procedure eliminates reactions which has substrates side empty
	 */
	public void cleanReactionsWithNoSubstrate() {
		List<Reaction> reactionsToRemove = new ArrayList<Reaction>();
		
		// for each reaction checks whether their products or substrates are empty
		for(Reaction r: getReactions().values()) {
			if(r.getSubstrates().isEmpty()){
				reactionsToRemove.add(r);
			}
		}

		for(Reaction rToDel: reactionsToRemove) {
			System.out.println("Without subs " + rToDel);
			removeReaction(rToDel);
		}
	}	

	/*
	 * This procedure eliminates repeated reactions 
	 */


	public void cleanRepeatedReactions() {
		List<Reaction> reactionsToRemove = new ArrayList<Reaction>();
		List<Reaction> reactions = new ArrayList<Reaction>(getReactions().values());

		// check repetitions and add them to be deleted
		for(int i = 0; i < reactions.size(); i++) {
			if(reactions.get(i).getId().endsWith("_REV")){
				continue;
			}

			// if the i-th reaction is already marked for removal, ignore it.
			if( reactionsToRemove.contains(reactions.get(i)) ) {
				continue;
			}



			// search repeated reactions and mark them
			for(int j = 0; j < reactions.size(); j++) {
				if(reactions.get(j).getId().endsWith("_REV") || i == j){
					continue;
				}

				if( reactions.get(i).haveSameSubstratesAndProducts(reactions.get(j)) ) {
					System.out.println("duplicate: " + reactions.get(i) + " and " + reactions.get(j));
					reactionsToRemove.add(reactions.get(j));
				}
			}
		}
		
		// delete the repeated reactions
		for(Reaction rToDel: reactionsToRemove) {
			//if(InputParameters.verbose){
			System.out.println("Removing repeated reaction "+rToDel);
			//}
			removeReactionAndItsReverse(rToDel);
		}
	}


	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public void cleanOrphanCompounds() {
		List<Compound> compoundsToRemove = new ArrayList<Compound>();

		// for each compound, check if it is not involved in any reaction
		for(Compound c: getCompounds().values()) {
			if( c.getProducedBy().isEmpty() && c.getSubstrateOf().isEmpty() ) {
				compoundsToRemove.add(c);
			}
		}

		for(Compound cToDel: compoundsToRemove) {
			removeCompound(cToDel, false);
		}
	}

	/*
	 * Initializes side compounds set with the product set for all reactions
	 * 
	 */
	public void initSideCompounds() {
		for(Reaction r: getReactions().values()) {
			r.initSideCompounds();
		}
	}

	public void removeCompoundsFromReactionsIfPresentAsSubstrateAndProduct() {
		for(Reaction r: getReactions().values()) {
			List<Compound> catalystCompounds = new ArrayList<Compound>();
			for(Compound s: r.getSubstrates().values()) {
				for(Compound p: r.getProduces().values() ) {
					if( p.getId().equals(s.getId()) ) {
						catalystCompounds.add(s);
						break;
					}
				}
			}

			// for each catalyst compound, remove it from both sides of the reaction
			for(Compound c: catalystCompounds) {
				r.removeSubstrate(c);
				r.removeProduct(c);
			}
		}
	}

	public String toString() {

		String str = new String();

		for(Reaction reaction : reactions.values()) {
			str = str.concat(reaction+"\n");
		}

		return str;
	}

	public HashMap<String, Compound> getCompounds() {
		return compounds;
	}

	public void setCompounds(HashMap<String, Compound> compounds) {
		this.compounds = compounds;
	}

	public int getNumPrecursors() {
		int numPrecursors = 0;
		for(Compound c: getCompounds().values()) {
			if( c.isPrecursor() )
				numPrecursors++;
		}
		return numPrecursors;
	}

	public int getNbInternalCompounds(){
		int nbIntComp = 0;
		for(Compound c: getCompounds().values()) {
			if( ! c.isPrecursor() )
				nbIntComp++;
		}
		return nbIntComp;
	}

	public HashMap<String, Reaction> getReactions() {
		return reactions;
	}

	public void setReactions(HashMap<String, Reaction> reactions) {
		this.reactions = reactions;
	}

	public void clearAllCompoundsFlags() {
		for(Compound c: compounds.values()) {
			c.clearFlag();
		}
	}

	public void clearAllReactionFlags() {
		for(Reaction r: reactions.values()) {
			r.clearFlag();
		}
	}

	public void removeUnBalancedReactions(){
		Iterator<Map.Entry<String, Reaction>> iterReac = this.getReactions().entrySet().iterator();
		while(iterReac.hasNext()){
			Reaction r = iterReac.next().getValue();
			/*HashMap<String, Integer> substratesAtoms = new HashMap<String, Integer>();
			HashMap<String, Integer> productsAtoms = new HashMap<String, Integer>();
			boolean substrateWithoutFormula = false;
			boolean productWithoutFormula = false;
			// get all atoms of substrates
			for(Compound subs : r.getSubstrates().values()){
				if(subs.getAtoms().keySet().size() == 0){
					substrateWithoutFormula = true;
					System.out.println(subs + " without formula");
				}
				for(String atomName : subs.getAtoms().keySet()){
					if(! substratesAtoms.containsKey(atomName)){
						substratesAtoms.put(atomName, (int) (r.getSubstratesStoich().get(subs.getId()) * subs.getAtoms().get(atomName)));
					}
					else{
						substratesAtoms.put(atomName, (int) (substratesAtoms.get(atomName) + r.getSubstratesStoich().get(subs.getId()) * subs.getAtoms().get(atomName)));
					}
				}
			}

			// get all atoms of products
			for(Compound prod : r.getProduces().values()){
				if(prod.getAtoms().keySet().size() == 0){
					productWithoutFormula = true;
					System.out.println(prod + " without formula");
				}
				for(String atomName : prod.getAtoms().keySet()){
					if(! productsAtoms.containsKey(atomName)){
						productsAtoms.put(atomName, (int) (r.getProductsStoich().get(prod.getId()) * prod.getAtoms().get(atomName)));
					}
					else{
						productsAtoms.put(atomName, (int) (productsAtoms.get(atomName) + r.getProductsStoich().get(prod.getId()) * prod.getAtoms().get(atomName)));
					}
				}
			}

			System.out.println(r + "\nSubstrates atoms sum: " + substratesAtoms);
			System.out.println("Products atoms sum: " + productsAtoms);
			if(substrateWithoutFormula || productWithoutFormula){
				System.out.println("Accept reaction because of compounds without formula");
				continue;
			}
			if(substratesAtoms.keySet().size() != productsAtoms.keySet().size()){
				System.out.println("Unbalanced in the original network ");
				iterReac.remove();
				removeReaction(r);
				continue;
			}
			for(String atomName : substratesAtoms.keySet()){
				if(! productsAtoms.containsKey(atomName) || substratesAtoms.get(atomName).intValue() != productsAtoms.get(atomName).intValue()){
					System.out.println("Unbalanced in the original network because of " + atomName);
					iterReac.remove();
					removeReaction(r);
					break;
				}
			}
			 */
			int[] nbCAtoms = new int[2]; 
			nbCAtoms[0] = nbCAtoms[1] = 0;
			boolean hasCompoundWithoutSmile = false;
			for(String subs : r.getSubstrates().keySet()){
				if(this.compounds.get(subs).getAtoms().keySet().size() == 0){
					hasCompoundWithoutSmile = true;
				}
				nbCAtoms[0] += r.getSubstrates().get(subs).getNumberOfAtom("C") * r.getSubstratesStoich().get(subs);
			}
			for(String prod : r.getProduces().keySet()){
				if(this.compounds.get(prod).getAtoms().keySet().size() == 0){
					hasCompoundWithoutSmile = true;
				}
				nbCAtoms[1] += r.getProduces().get(prod).getNumberOfAtom("C") * r.getProductsStoich().get(prod);
			}
			if(hasCompoundWithoutSmile && r.getSubstrates().keySet().size() > 0 && r.getProduces().keySet().size() > 0){
				if(InputParameters.verbose){
					System.out.println("contain comp without smiles: " + r);
				}
				continue;
			}
			if(nbCAtoms[0] != nbCAtoms[1]){ // not balanced
				if(InputParameters.verbose){
					System.out.println("Unbalanced in the original network " + r);
				}
				iterReac.remove();
				removeReaction(r);
			}

		}
	}

	public void removePairedCoFactors(boolean removeCoA) {
		int nbNewReactions = 0;
		List<String> newCofReactions = new ArrayList<String>();
		// check for every reaction if we can remove paired cofactors
		Iterator<Map.Entry<String, Reaction>> iterReac = this.getReactions().entrySet().iterator();
		while(iterReac.hasNext()){
			Reaction r = iterReac.next().getValue();
			if(r.getCofactors().size() == 0 || r.getId().endsWith("_REV")){ // no cofactors
				continue;
			}
			//System.out.print("reaction " + r.getId());
			// check if we have still at least one substrate and one product
			// check if the reaction is still balanced (C-atoms)
			// first two entries -> nb of C-Atoms without co-factors
			// last two entries -> nb of C-Atoms of all substrates and products
			int[] nbCAtoms = new int[4]; 
			nbCAtoms[0] = nbCAtoms[1] = nbCAtoms[2] = nbCAtoms[3] = 0;
			boolean hasOtherSubstratesThanCoFactors = false;
			boolean hasOtherProductsThanCoFactors = false;
			boolean hasCompoundWithoutSmile = false;
			for(Compound subs : r.getSubstrates().values()){
				// substrate without chem. formula
				if(subs.getAtoms().keySet().size() == 0){
					hasCompoundWithoutSmile = true;
				}
				if(! r.getCofactors().contains(subs)){
					nbCAtoms[0] += r.getSubstrates().get(subs.getId()).getNumberOfAtom("C") * r.getSubstratesStoich().get(subs.getId());
					hasOtherSubstratesThanCoFactors = true;

				}
				nbCAtoms[2] += r.getSubstrates().get(subs.getId()).getNumberOfAtom("C") * r.getSubstratesStoich().get(subs.getId());
			}
			for(Compound prod : r.getProduces().values()){
				// product without chem. formula
				if(prod.getAtoms().keySet().size() == 0){
					hasCompoundWithoutSmile = true;
				}
				if(! r.getCofactors().contains(prod)){
					nbCAtoms[1] += r.getProduces().get(prod.getId()).getNumberOfAtom("C") * r.getProductsStoich().get(prod.getId());
					hasOtherProductsThanCoFactors = true;
				}
				nbCAtoms[3] += r.getProduces().get(prod.getId()).getNumberOfAtom("C") * r.getProductsStoich().get(prod.getId());
			}
			if(hasCompoundWithoutSmile){
				if(InputParameters.verbose){
					System.out.println("cof contain compounds without smiles: " + r);
				}
				continue;
			}
			// the reaction is unbalanced in the original network -> remove it
			if(nbCAtoms[2] != nbCAtoms[3]){
				if(InputParameters.verbose){
					System.out.println("Unbalanced in the original network " + r);
				}

				iterReac.remove();
				removeReaction(r);
				continue;
			}


			// reaction is balanced after removal of co-factors
			// split co-factors and build a reaction between them
			if(hasOtherSubstratesThanCoFactors && hasOtherProductsThanCoFactors && nbCAtoms[0] == nbCAtoms[1]){

				//if(nbCAtoms[0] == nbCAtoms[1]){
				//System.out.println("...balanced");
				if(InputParameters.verbose){
					System.out.println("with cofactors " + r);
				}

				// remove co-factors 
				List<String> newReactions = new ArrayList<String>();
				newReactions.addAll(buildCofactorReactions(r, nbNewReactions, true));
				newCofReactions.addAll(newReactions);
				nbNewReactions += newReactions.size();
				// remove co-factors from original reaction
				for(Compound c : r.getCofactors()){
					if(r.getProduces().containsKey(c.getId())){
						r.removeProduct(c);
					}
					else if(r.getSubstrates().containsKey(c.getId())){
						r.removeSubstrate(c);
					}
				}
				if(InputParameters.verbose){
					System.out.println("Reaction without cofactors " + r);
				}

			}
			// reaction contains only co-factors -> try to split the co-factor pairs
			// e.g.: ATP + NAD -> ADP + NADP becomes ATP -> ADP and NAD -> NADP
			else if(! hasOtherSubstratesThanCoFactors && ! hasOtherProductsThanCoFactors && nbCAtoms[0] == nbCAtoms[1]){
				if(InputParameters.verbose){
					System.out.println("Reaction with only cofactors " + r);
				}
				// try to find the right split
				List<String> newReactions = new ArrayList<String>();
				newReactions.addAll(buildCofactorReactions(r, nbNewReactions, true));
				newCofReactions.addAll(newReactions);
				nbNewReactions += newReactions.size();

				// delete old reaction that contains only (paired) cofactors
				iterReac.remove();
				removeReaction(r);
			}
			// not balanced if we remove all co-factors of the reaction,
			// but can we remove a part of the co-factors?
			else{
				if(InputParameters.verbose){
					System.out.println("reaction " + r + "...not balanced if we try to remove " + r.getCofactors() );
				}


				if(r.getSubstrates().keySet().size() > 2 || r.getProduces().keySet().size() > 2){
					if(InputParameters.verbose){
						for(Compound c : r.getSubstrates().values()){
							System.out.println("Subs " + c.getId() + " " + c.getNumberOfAtom("C") + " x " + r.getSubstrateStochiometricValue(c));
						}
						for(Compound c : r.getProduces().values()){
							System.out.println("Prod " + c.getId() + " " + c.getNumberOfAtom("C") + " x " + r.getProductStochiometricValue(c));
						}
						System.out.println("Try to remove part of the co-factors");
					}
					List<String> newReactions = new ArrayList<String>();
					newReactions.addAll(buildCofactorReactions(r, nbNewReactions, false));
					if(newReactions.size() > 0){
						Iterator<String> iterStr = newReactions.iterator();
						while(iterStr.hasNext()){
							String reacString = iterStr.next();
							String copy = new String(reacString);
							copy = copy.concat(r.getId());
							newCofReactions.add(copy);
						}
					}
					//newCofReactions.addAll(newReactions);
					nbNewReactions += newReactions.size();

					/*System.out.println("New reactions: " + newReactions);
					if(newReactions.size() > 0){
						for(Reaction reac : newReactions){
							for(Compound c : reac.getSubstrates().values()){
								r.removeSubstrate(c);
							}
							for(Compound c : reac.getProduces().values()){
								r.removeProduct(c);
							}
						}
						if(r.getSubstrates().keySet().size() == 0 || r.getProduces().keySet().size() == 0){
							System.out.println("Delete empty reaction " + r);
							iterReac.remove();
							removeReaction(r);
						}
						else{
							System.out.println("Without cofactors " + r);
						}
					}*/
				}

				/*if(removeCoA){
					removeCoA(r);
				}*/
			}

			if(InputParameters.verbose){
				System.out.println("----------------------");
			}

		}


		// add all co-factor reactions to a hash to filter duplicates
		Map<String, String> groupCofReations = new HashMap<String, String>();
		for(String reacString : newCofReactions){
			if(InputParameters.verbose){
				System.out.println("New cof: " + reacString);
			}
			// 0.. id; 1..reversibility; 2..substrates; 3..products; 4.. (optional) remove compounds from this reaction
			String[] parts = reacString.split("###");
			if(InputParameters.verbose){
				System.out.println(Arrays.toString(parts));
			}
			String key = parts[1] + "###" + parts[2] + "###" + parts[3];
			// add the reaction where we have to split off the cof
			if(! groupCofReations.containsKey(key) || parts.length == 5){ 
				String reac = "";
				if(parts.length == 5){
					if(groupCofReations.containsKey(key)){
						reac = groupCofReations.get(key);
					}
					reac = reac.concat("##" + parts[4]);
				}
				groupCofReations.put(key, reac);
			}
		}
		if(InputParameters.verbose){
			for(String s : groupCofReations.keySet()){
				System.out.println("map " + s + " -> " + groupCofReations.get(s));
			}
		}
		// build co-factor reactions from strings
		int idx = 0;
		for(String reacString : groupCofReations.keySet()){
			if(InputParameters.verbose){
				System.out.println("treat " + reacString + " -> " + groupCofReations.get(reacString));
			}
			idx++;
			String[] parts = reacString.split("###");
			boolean reversible = false;
			if(parts[0].equals("+")){
				reversible = true;
			}
			Reaction r1 = addNewReaction("COF__45__" + idx, "COF__45__" + idx, reversible);
			Reaction r2 = null;
			if(reversible){
				r2 = addNewReaction("COF__45__" + idx + "_REV", "COF__45__" + idx + "_REV", reversible);
				r1.setReverse(r2);
			}


			// add substrates
			//Map<String, Double> substrates = new HashMap<String, Double>();
			String[] subs = parts[1].split("##");
			for(String s : subs){
				String[] items = s.split("#");
				if(InputParameters.verbose){
					System.out.println("Subs: " + Arrays.toString(items));
				}
				for(int i = 0; i < items.length; i += 2){
					//substrates.put(items[i], Double.valueOf(items[i+1]));
					r1.addSubstrate(this.compounds.get(items[i]), Double.valueOf(items[i+1]));
					if(reversible){
						r2.addProduct(this.compounds.get(items[i]), Double.valueOf(items[i+1]));
					}
				}
			}
			//System.out.println("Subs:" + substrates);

			// add products
			//Map<String, Double> products = new HashMap<String, Double>();
			String[] prod = parts[2].split("##");
			for(String s : prod){
				String[] items = s.split("#");
				if(InputParameters.verbose){
					System.out.println("Prod: " + Arrays.toString(items));
				}
				for(int i = 0; i < items.length; i += 2){
					//products.put(items[i], Double.valueOf(items[i+1]));
					r1.addProduct(this.compounds.get(items[i]), Double.valueOf(items[i+1]));
					if(reversible){
						r2.addSubstrate(this.compounds.get(items[i]), Double.valueOf(items[i+1]));
					}
				}
			}

			// normalise stoichiometry
			double stoichiometryFactor = 1.0;
			List<Double> stoichiometricEntries = new ArrayList<Double>();
			for(Compound c : r1.getSubstrates().values()){
				stoichiometricEntries.add(r1.getSubstrateStochiometricValue(c));
			}
			for(Compound c : r1.getProduces().values()){
				stoichiometricEntries.add(r1.getProductStochiometricValue(c));
			}
			Double minElem = Collections.min(stoichiometricEntries);
			if(minElem > 1.0){
				if(InputParameters.verbose){
					System.out.println("Higher stoich " + minElem);
				}
				boolean stillIntValue = true;
				for(Double d : stoichiometricEntries){
					if(((d/minElem) % 1) != 0){
						stillIntValue = false;
						break;
					}
				}
				if(stillIntValue == true){
					if(InputParameters.verbose){
						System.out.println("Divide all entries by " + minElem);
					}
					stoichiometryFactor = minElem;
				}
			}

			Iterator<Compound> iterComp = r1.getSubstrates().values().iterator();
			while(iterComp.hasNext()){
				Compound c = iterComp.next();
				r1.addSubstrate(c, r1.getSubstrateStochiometricValue(c) / stoichiometryFactor);
				if(r1.isReversible()){
					r2.addProduct(c, r1.getSubstrateStochiometricValue(c) / stoichiometryFactor);
				}
			}


			iterComp = r1.getProduces().values().iterator();
			while(iterComp.hasNext()){
				Compound c = iterComp.next();
				r1.addProduct(c, r1.getProductStochiometricValue(c) / stoichiometryFactor);
				if(r1.isReversible()){
					r2.addSubstrate(c, r1.getProductStochiometricValue(c) / stoichiometryFactor);
				}
			}

			// remove cofactors from original reaction
			if(! groupCofReations.get(reacString).equals("")){
				String[] removeCofFromReac = groupCofReations.get(reacString).split("##");
				for(String removePartOfReaction : removeCofFromReac){
					if(removePartOfReaction.equals("")){
						continue;
					}

					Reaction r = this.reactions.get(removePartOfReaction);
					if(InputParameters.verbose){
						System.out.println("Remove cofactors from this reaction: " + removePartOfReaction);
						System.out.println("r " + r);
					}
					for(Compound c : r1.getSubstrates().values()){
						r.removeSubstrate(c);
						if(r.isReversible()){
							r.getReverseReaction().removeProduct(c);
						}
					}
					for(Compound c : r1.getProduces().values()){
						r.removeProduct(c);
						if(r.isReversible()){
							r.getReverseReaction().removeSubstrate(c);
						}
					}
					if(InputParameters.verbose){
						System.out.println("Without cofactors " + r);
					}
				}
			}
		}

		iterReac = this.getReactions().entrySet().iterator();
		while(iterReac.hasNext()){
			Reaction r = iterReac.next().getValue();
			if(removeCoA){
				removeCoA(r);
			}
		}
	}

	private void removeCoA(Reaction r){
		if(r.getSubstrates().containsKey("CO__45__A") && r.getProduces().containsKey("ACETYL__45__COA") &&
				r.getSubstrates().keySet().size() > 1 && r.getProduces().keySet().size() > 1){
			//System.out.println("Remove coA from " + r);
			Iterator<Compound> iterComp = r.getSubstrates().values().iterator();
			while(iterComp.hasNext()){
				Compound c = iterComp.next();
				if(c.getId().equals("CO__45__A")){
					iterComp.remove();
					r.removeSubstrate(c);
				}
			}
			iterComp = r.getProduces().values().iterator();
			while(iterComp.hasNext()){
				Compound c = iterComp.next();
				if(c.getId().equals("ACETYL__45__COA")){
					iterComp.remove();
					r.removeProduct(c);
				}
			}
		}
		else if(r.getSubstrates().containsKey("ACETYL__45__COA") && r.getProduces().containsKey("CO__45__A") &&
				r.getSubstrates().keySet().size() > 1 && r.getProduces().keySet().size() > 1){
			//System.out.println("Remove coA from " + r);
			Iterator<Compound> iterComp = r.getSubstrates().values().iterator();
			while(iterComp.hasNext()){
				Compound c = iterComp.next();
				if(c.getId().equals("ACETYL__45__COA")){
					iterComp.remove();
					r.removeSubstrate(c);
				}
			}
			iterComp = r.getProduces().values().iterator();
			while(iterComp.hasNext()){
				Compound c = iterComp.next();
				if(c.getId().equals("CO__45__A")){
					iterComp.remove();
					r.removeProduct(c);
				}
			}
		}
		else if((r.getSubstrates().containsKey("CO__45__A") && r.getProduces().containsKey("ACETYL__45__COA")) ||
				(r.getSubstrates().containsKey("ACETYL__45__COA") && r.getProduces().containsKey("CO__45__A"))){
			System.out.println("can not delete coA from " + r);
		}
	}

	// return a list of strings, where each string stays for a reaction that will be created later
	// a string is as follow: id###[+|-]###subs_1#stoich_1##...##subs_n#stoich_n###prod_1#stoich_1##...##prod_m#stoich_m
	private List<String> buildCofactorReactions(Reaction r, int counter, boolean onlyCofactors) {
		List<String> returnList = new ArrayList<String>();
		List<Compound> substrateList = new ArrayList<Compound>();
		List<Compound> productList = new ArrayList<Compound>();
		List<Compound> cmpList = new ArrayList<Compound>();
		// build list of substrates and products of the reactions
		// consider only co-factors or all compounds
		if(onlyCofactors == true){
			for(Compound c : r.getCofactors()){
				cmpList.add(c);
			}
		}
		else{
			for(Compound c : r.getSubstrates().values()){
				cmpList.add(c);
			}
			for(Compound c : r.getProduces().values()){
				cmpList.add(c);
			}
		}
		for(Compound c : cmpList){
			if(r.getSubstrates().containsKey(c.getId())){
				substrateList.add(c);	
			}
			else{
				productList.add(c);
			}
		}
		// build the power set (without the empty set) of the compound list 
		// associated with the sum of 'C' atoms to each combination
		LinkedHashMap<List<Compound>, Double> substrateSum = sumAtomsOfAllCombinations(r, substrateList);
		LinkedHashMap<List<Compound>, Double> productSum = sumAtomsOfAllCombinations(r, productList);

		// filter known co-factor pairs
		// here we could also take a user-provided list...
		List<List<Compound>> substratesOfValidatedNewReactions = new ArrayList<List<Compound>>();
		List<List<Compound>> productsOfValidatedNewReactions = new ArrayList<List<Compound>>();
		Iterator<List<Compound>> iter1 = substrateSum.keySet().iterator();
		while(iter1.hasNext()){
			List<Compound> l1 = iter1.next();
			Iterator<List<Compound>> iter2 = productSum.keySet().iterator();
			boolean foundPair = false;
			while(iter2.hasNext()){
				List<Compound> l2 = iter2.next();
				if(l1.size() == 1 && l2.size() == 1 && (
						(l1.get(0).getId().equals("NAD") && l2.get(0).getId().equals("NADH")) ||
						(l1.get(0).getId().equals("NADH") && l2.get(0).getId().equals("NAD")) 		
						)){
					// create reaction
					// id == name
					String newReac = new String("COF__45__" + counter);
					// reversibility
					if(r.isReversible()){
						newReac = newReac.concat("###+");
					}
					else{
						newReac = newReac.concat("###-");
					}
					// substrates
					newReac = newReac.concat("###" + l1.get(0).getId() + "#" + r.getSubstrateStochiometricValue(l1.get(0)));
					//products
					newReac = newReac.concat("###" + l2.get(0).getId() + "#" + r.getProductStochiometricValue(l2.get(0)) + "###");
					if(InputParameters.verbose){
						System.out.println("Build reaction " + newReac);
					}
					returnList.add(newReac);
					// reverse reaction ?




					counter++;
					substratesOfValidatedNewReactions.add(l1);
					productsOfValidatedNewReactions.add(l2);

					// delete items from lists
					iter2.remove();
					foundPair = true;
					break;
				}
			}
			if(foundPair == true){
				iter1.remove();
			}
		}

		if(InputParameters.verbose){
			System.out.println("HashMap subs: " + substrateSum);
			System.out.println("HashMap prod: " + productSum);
		}

		// find minimal combinations of substates and products that have the same number of 'C' atoms
		// if we can build A -> C, B -> D, and A+B -> C+D we build only the first two reactions
		for(List<Compound> subsComposition : substrateSum.keySet()){
			// exclude supersets of substrates of created reaction 
			boolean alreadyASubstrate = false;
			for(List<Compound> item : substratesOfValidatedNewReactions){
				if(subsComposition.containsAll(item)){
					alreadyASubstrate = true;
					break;
				}
			}
			if(alreadyASubstrate == true){
				continue;
			}

			List<Compound> prodCompositionOfNewReaction = new ArrayList<Compound>();
			int nbMatches = 0;
			//System.out.println("subs " + subsComposition + " -> " + substrateSum.get(subsComposition));
			for(List<Compound> prodComposition : productSum.keySet()){
				// check if products are superset of former solution
				boolean alreadyAProduct = false;
				for(List<Compound> item : productsOfValidatedNewReactions){
					if(prodComposition.containsAll(item)){
						alreadyAProduct = true;
						break;
					}
				}
				if(alreadyAProduct == true){
					continue;
				}
				//System.out.println("prod " + prodComposition + " -> " + productSum.get(prodComposition));
				if(substrateSum.get(subsComposition).equals(productSum.get(prodComposition))){
					nbMatches++;
					prodCompositionOfNewReaction = prodComposition;
					//	System.out.println("Build reaction between " + subsComposition + " and " + prodComposition);
				}
			}
			// build new reaction
			if(nbMatches == 1){
				substratesOfValidatedNewReactions.add(substrateList);
				String newReac = new String("COF__45__" + counter);
				if(r.isReversible()){
					newReac = newReac.concat("###+###");
				}
				else{
					newReac = newReac.concat("###-###");
				}

				//Reaction cofactorReac = addNewReaction("COF__45__" + counter, "COF__45__" + counter, r.isReversible());

				// add substrates and products
				for(Compound c : r.getSubstrates().values()){
					if(subsComposition.contains(c)){
						//cofactorReac.addSubstrate(c, r.getSubstrateStochiometricValue(c));
						newReac = newReac.concat(c.getId() + "#" + r.getSubstrateStochiometricValue(c) + "##");
					}
				}
				newReac = newReac.concat("#");
				for(Compound c : r.getProduces().values()){
					if(prodCompositionOfNewReaction.contains(c)){
						//cofactorReac.addProduct(c, r.getProductStochiometricValue(c));
						newReac = newReac.concat(c.getId() + "#" + r.getProductStochiometricValue(c) + "##");
					}
				}
				newReac = newReac.concat("#");
				//System.out.println("Build reaction " + cofactorReac);
				if(InputParameters.verbose){
					System.out.println("Build reaction " + newReac);
				}
				returnList.add(newReac);


				counter++;
				//	System.out.println("Build new reaction " + cofactorReac);
			}
		}

		return returnList;
	}

	private LinkedHashMap<List<Compound>, Double> sumAtomsOfAllCombinations(Reaction r, List<Compound> cmpList){
		if(r.getId().equals("ACETALD__45__DEHYDROG__45__RXN")){
			System.out.println("Check");
		}
		LinkedHashMap<List<Compound>, Double> sumMap = new LinkedHashMap<List<Compound>, Double>();
		int binSize = (int)Math.pow(2, cmpList.size());
		for(Integer i= 1;i< binSize;++i){
			String bin= Integer.toString(i, 2); //convert to binary
			while(bin.length() < cmpList.size())bin = "0" + bin; //pad with 0's
			List<Compound> thisComb = new ArrayList<Compound>(); //place to put one combination
			Double nbCAtoms = 0.0;
			for(int j= 0;j< cmpList.size();++j){
				if(bin.charAt(j) == '1'){
					thisComb.add(cmpList.get(j));
					if(r.substrates.containsKey(cmpList.get(j).getId())){
						nbCAtoms += r.getSubstrateStochiometricValue(cmpList.get(j)) * cmpList.get(j).getNumberOfAtom("C");
					}
					else{
						nbCAtoms += r.getProductStochiometricValue(cmpList.get(j)) * cmpList.get(j).getNumberOfAtom("C");
					}
				}
			}
			if(thisComb.size() == 0){
				continue;
			}
			Collections.sort(thisComb); //sort it for easy checking
			sumMap.put(thisComb, nbCAtoms); //put this set in the answer list
		}

		return sumMap;
	}


	public boolean checkNetwork(List<Compound> targets, List<Compound> compoundList){
		System.out.println("Remove reactions that allows a compound to be produced from any source");
		// check if we can produce each non-source compound from the union of all sources
		List<Compound> allSources = new ArrayList<Compound>();
		for(Compound c : this.compounds.values()){
			if(c.isPrecursor()){
				allSources.add(c);
			}
		}


		List<Reaction> networkReactions = new ArrayList<Reaction>();
		for (Reaction reac : this.reactions.values()) {
			networkReactions.add(reac);
		}
		// do the checks
		int nbInconsistance = 0;
		// delete reactions that produce compounds out of nothing
		Map<Reaction, List<Compound>> delReactions = new HashMap<Reaction, List<Compound>>();

		List<Compound> noSourceList = new ArrayList<Compound>();
		for(Compound c : this.compounds.values()){
			if(! c.isPrecursor()){
				noSourceList.add(c);
			}
		}

		for(Reaction r : delReactions.keySet()){
			if(InputParameters.verbose){
				System.out.println("Remove reaction " + r.getId() + " because it allows the production from any sorce of these compounds -> " + delReactions.get(r));
			}
			removeReaction(r);
		}



		//System.exit(0);
		if(nbInconsistance > 0){
			return false;
		}
		else{
			return true;
		}
	}


	public void parseSbmlFormat(String sbmlFile) {
		// TODO Auto-generated method stub
		SAXBuilder sxb = new SAXBuilder();
		try {
			Document document = sxb.build(new File(sbmlFile));
			Element root = document.getRootElement();
			addCompounds(root);
			addReactions(root);
		} catch (JDOMException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void addCompounds(Element root){
		for(Element elem : root.getChildren()){
			if(elem.getName().equals("model")){
				for(Element elem2 : elem.getChildren()){
					if(elem2.getName().equals("listOfSpecies")){
						for(Element species : elem2.getChildren()){
							if(species.getName().equals("species")){
								Compound c = addCompound(species.getAttributeValue("id"), species.getAttributeValue("name"), 
										species.getAttributeValue("compartment"));
								String boundaryString = species.getAttributeValue("boundaryCondition");
								if (boundaryString != null){
									System.out.print("\"" + c.getId() + "\"" + ",");
									c.setBoundary(boundaryString.toLowerCase().equals("true") ? true : false);
								}
								// get atoms and number of atoms 
								// build hash map AtomName => #atoms
								HashMap<String, Integer> atomMap = new HashMap<String, Integer>();
								String formula = species.getAttributeValue("formula");
								if(formula != null && ! formula.equals("")){
									String[] atoms = formula.split("\\d+"); // split on numbers
									String[] nbAtoms = formula.split("\\D+"); // split on letters
									for(int i = 0; i < atoms.length; i++){
										atomMap.put(atoms[i], Integer.parseInt(nbAtoms[i+1])); // nbAtoms[i+1], because nbAtoms[0] is empty -> string before formula
									}
								}
								c.setAtoms(atomMap);
							}
						}
					}
				}
			}
		}	
		System.out.println("");
	}

	private void addReactions(Element root){
		for(Element elem : root.getChildren()){
			if(elem.getName().equals("model")){
				for(Element elem2 : elem.getChildren()){
					if(elem2.getName().equals("listOfReactions")){
						for(Element reaction : elem2.getChildren()){
							if(reaction.getName().equals("reaction")){
								Reaction r1 = addNewReaction(reaction.getAttributeValue("id"), 
										reaction.getAttributeValue("name"), Boolean.valueOf(reaction.getAttributeValue("reversible")));
								Reaction r2 = null;
								if( r1.isReversible() ) {
									r2 = addNewReaction(r1.getId()+"_REV", r1.getName()+"_REV", r1.isReversible());
									r1.setReverse(r2);
								}

								// add the substrates/products/cofactors of the reaction
								for(Element elem3 : reaction.getChildren()){
									if(elem3.getName().equals("listOfReactants")){
										for(Element speciesRef : elem3.getChildren()){
											String cId = speciesRef.getAttributeValue("species");
											double stoich = 1.0;
											try{
												stoich = Double.valueOf(speciesRef.getAttributeValue("stoichiometry"));
											}
											catch (Exception e){
												System.err.println("Stoichiometry value not present, assuming 1.0 .");
											}
											Compound c = compounds.get( cId );
											if( c != null)
											{
												r1.addSubstrate(c, stoich);
												if( r1.isReversible() )
													r2.addProduct(c, stoich);
											}
										}
									}
									else if(elem3.getName().equals("listOfProducts")){
										for(Element speciesRef : elem3.getChildren()){
											String cId = speciesRef.getAttributeValue("species");
											double stoich = 1.0;
											try{
												stoich = Double.valueOf(speciesRef.getAttributeValue("stoichiometry"));
											}
											catch (Exception e){
												System.err.println("Stoichiometry value not present, assuming 1.0 .");
											}
											Compound c = compounds.get( cId );
											if( c != null)
											{
												r1.addProduct(c, stoich);
												if( r1.isReversible() )
													r2.addSubstrate(c, stoich);
											}
										}
									}
									else if(elem3.getName().equals("cofactors")){
										for(Element speciesRef : elem3.getChildren()){
											String cId = speciesRef.getAttributeValue("species");
											Compound c = compounds.get( cId );
											if( c != null){
												r1.addCofactor(c);
												if(r1.isReversible())
													r2.addCofactor(c);
											}
										}
									}
								}
								if(r1.getId().equals("R_GLCNtex")){
									System.out.println();
								}
							}
						}
					}
				}
			}
		}
	}


	public void markHighDegreeCompounds() {
		// TODO Auto-generated method stub
		HashMap<Integer, List<Compound>> incomingDeg = new HashMap<Integer, List<Compound>>();
		HashMap<Integer, List<Compound>> outgoingDeg = new HashMap<Integer, List<Compound>>();
		for(Compound c : this.getCompounds().values()){
			if(c.isPrecursor() || InputParameters.getTargetCompounds().contains(c)){
				continue;
			}
			int key = c.getProducedBy().size();
			List<Compound> l;
			if(incomingDeg.containsKey(key)){
				l = incomingDeg.get(key);
			}
			else{
				l = new ArrayList<Compound>();
			}
			l.add(c);
			incomingDeg.put(key, l);

			key = c.getSubstrateOf().size();
			if(outgoingDeg.containsKey(key)){
				l = outgoingDeg.get(key);
			}
			else{
				l = new ArrayList<Compound>();
			}
			l.add(c);
			outgoingDeg.put(key, l);
		}

		Set<Compound> highDegreeCompounds = new HashSet<Compound>();
		for(Integer degree : incomingDeg.keySet()){
			if(degree >= 10){
				highDegreeCompounds.addAll(incomingDeg.get(degree));
				System.out.println("Degree: " + degree + " -> " + incomingDeg.get(degree));
			}
		}
		System.out.println("---------");
		for(Integer degree : outgoingDeg.keySet()){
			if(degree >= 10){
				highDegreeCompounds.addAll(outgoingDeg.get(degree));
				System.out.println("Degree: " + degree + " -> " + outgoingDeg.get(degree));
			}
		}

		// mark all high degree compounds
		for(Compound c : highDegreeCompounds){
			c.setHighDegree(true);
		}
	}


	public void removeCO2() {
		// TODO Auto-generated method stub
		for(Reaction r : this.reactions.values()){
			if(InputParameters.verbose){
				System.out.println("Remove co2 from " + r);
			}
			Iterator<Compound> iterSubs = r.getSubstrates().values().iterator();
			Compound co2Compound = null;
			while(iterSubs.hasNext()){
				Compound co2 = iterSubs.next();
				if(co2.getId().equals("CARBON__45__DIOXIDE")){
					co2Compound = co2;
					//iterSubs.remove();
					//r.removeSubstrate(co2);

				}
			}
			if(co2Compound != null){
				r.removeSubstrate(co2Compound);
			}

			co2Compound = null;
			Iterator<Compound> iterProd = r.getProduces().values().iterator();
			while(iterProd.hasNext()){
				Compound co2 = iterProd.next();
				if(co2.getId().equals("CARBON__45__DIOXIDE")){
					co2Compound = co2;
					//iterProd.remove();
					//r.removeProduct(co2);

				}
			}
			if(co2Compound != null){
				r.removeProduct(co2Compound);
			}
			if(r.getId().equals("COF__45__107")){
				System.out.println();
			}
		}
		Iterator<String> iterComp = this.compounds.keySet().iterator();
		while(iterComp.hasNext()){
			String cId = iterComp.next();
			if(cId.equals("CARBON__45__DIOXIDE")){
				iterComp.remove();
				this.removeCompound(this.compounds.get(cId), false);

				break;
			}
		}
	}

	public void removeElectronAcceptors() {
		// TODO Auto-generated method stub
		for(Reaction r : this.reactions.values()){
			Iterator<Compound> iterSubs = r.getSubstrates().values().iterator();
			while(iterSubs.hasNext()){
				Compound ea = iterSubs.next();
				if(ea.getId().equals("Acceptor") || ea.getId().equals("Donor__45__H2")){
					iterSubs.remove();
					r.removeSubstrate(ea);

				}
			}

			Iterator<Compound> iterProd = r.getProduces().values().iterator();
			while(iterProd.hasNext()){
				Compound ea = iterProd.next();
				if(ea.getId().equals("Acceptor") || ea.getId().equals("Donor__45__H2")){
					iterProd.remove();
					r.removeProduct(ea);

				}
			}
		}
	}

	// if there is a reaction r2 that is subset of reaction r1
	// than can we remove all substrates and products from r1 that appears in r2
	public void removeSubReactions(){
		HashMap<Reaction, Set<Compound>> compoundsToRemovePerReaction = new HashMap<Reaction, Set<Compound>>();
		for(Reaction r1 : this.reactions.values()){
			for(Reaction r2 : this.reactions.values()){
				if(r1.equals(r2)){
					continue;
				}
				// if we have r1: A + B <-> C + D and r2: A -> C
				// we can not cut A,C from r2 
				if(! r2.isReversible() && r1.isReversible()){ 
					continue;
				}
				// reactions without substrates or products will be deleted in another step
				if(r1.getSubstrates().keySet().size() == 0 || r2.getSubstrates().keySet().size() == 0 ||
						r1.getProduces().values().size() == 0 || r2.getProduces().values().size() == 0){
					continue;
				}
				// r2 can not be a strict subset if it has at least as many substrates and products as r1
				if(r1.getSubstrates().keySet().size() <= r2.getSubstrates().keySet().size() ||
						r1.getProduces().keySet().size() <= r2.getProduces().keySet().size()){
					continue;
				}

				Set<Compound> compoundsToRemove = new HashSet<Compound>();
				boolean isSubset = true;
				for(Compound subs : r2.getSubstrates().values()){
					if(! r1.getSubstrates().containsKey(subs.getId()) || r1.getSubstrateStochiometricValue(subs) != r2.getSubstrateStochiometricValue(subs)){
						isSubset = false;
						break;
					}
					compoundsToRemove.add(subs);
				}

				if(isSubset == false){
					continue;
				}

				for(Compound prod : r2.getProduces().values()){
					if(! r1.getProduces().containsKey(prod.getId()) || r1.getProductStochiometricValue(prod) != r2.getProductStochiometricValue(prod)){
						isSubset = false;
						break;
					}
					compoundsToRemove.add(prod);
				}

				if(isSubset == false){
					continue;
				}

				// here we have subsets in sources and products
				// remove sources and products from r1
				if(compoundsToRemovePerReaction.containsKey(r1)){
					compoundsToRemove.addAll(compoundsToRemovePerReaction.get(r1));
				}
				compoundsToRemovePerReaction.put(r1, compoundsToRemove);
				if(InputParameters.verbose){
					System.out.println(r2 + " is subset of " + r1);
				}
			}
		}

		for(Reaction r : compoundsToRemovePerReaction.keySet()){
			if(InputParameters.verbose){
				System.out.println("Remove subsets " + compoundsToRemovePerReaction.get(r) + " from " + r);
			}
			// try to remove compound as product or substrate
			for(Compound c : compoundsToRemovePerReaction.get(r)){
				r.removeSubstrate(c);
				r.removeProduct(c);
			}
			if(InputParameters.verbose){
				System.out.println("New reaction: " + r);
			}
		}
	}


	public void checkForSubReactions() {
		// TODO Auto-generated method stub

		for(Reaction r1 : this.reactions.values()){
			for(Reaction r2 : this.reactions.values()){
				// same reaction or different number of products
				if(r1.equals(r2) || r1.getProduces().keySet().size() != r2.getProduces().keySet().size()){
					continue;
				}

				// check if they have the same products
				List<Compound> products = new ArrayList<Compound>();
				products.addAll(r1.getProduces().values());
				products.removeAll(r2.getProduces().values());
				if(products.size() != 0){
					continue;
				}

				if(r1.getSubstrates().keySet().size() >= r2.getSubstrates().keySet().size()){
					continue;
				}
				List<Compound> substrates = new ArrayList<Compound>();
				substrates.addAll(r1.getSubstrates().values());
				substrates.removeAll(r2.getSubstrates().values());
				if(substrates.size() == 0){
					System.out.println(r1 + " is a subset reaction of " + r2);
				}

			}
		}
	}


	public void countProducts() {
		// TODO Auto-generated method stub
		//init
		int maxOcc = 0;
		for(Reaction r : this.reactions.values()){
			if(maxOcc < r.getProduces().keySet().size()){
				maxOcc = r.getProduces().keySet().size();
			}
		}
		List<Integer> countNbProductOcc = new ArrayList<Integer>(maxOcc);
		System.out.println("Max number of products: " + maxOcc);
		for(int i = 0; i <= maxOcc; i++){
			countNbProductOcc.add(i, 0);
		}
		for(Reaction r : this.reactions.values()){
			countNbProductOcc.set(r.getProduces().keySet().size(), countNbProductOcc.get(r.getProduces().keySet().size()) + 1);
		}

		if(InputParameters.verbose){
			for(int i = 1; i <= maxOcc; i++){
				System.out.println("Nb of reactions with " + i + " product: " + countNbProductOcc.get(i) + " / " + this.reactions.keySet().size() + " -> " + (countNbProductOcc.get(i) / this.reactions.keySet().size()));
			}
		}
		//System.exit(0);
	}


	public void getUnproducableCompounds() {
		// TODO Auto-generated method stub

	}




}	




