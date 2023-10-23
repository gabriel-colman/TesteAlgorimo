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
package application;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import metabolicNetwork.Compound;
import metabolicNetwork.Reaction;

/**
 * @author ludo
 * This class represents a solution of the findPrecursors method. 
 * It contains the set of precursors and the compounds used as bootstrap in the auto fed cycles.
 * 
 *
 */
public class PrecursorSet implements Comparable<PrecursorSet>, java.io.Serializable{
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -8613887645675257673L;
	private List<Compound> precursors;
	private List<Compound> bootstraps;
	private Set<Reaction> reactions;
	private Set<Compound> cumulatedCompounds = new HashSet<Compound>();
	

	/**
	 * @return the cumullatedCompounds
	 */
	public Set<Compound> getCumulatedCompounds() {
		return cumulatedCompounds;
	}

	/**
	 * @param cumullatedCompounds the cumulatedCompounds to set
	 */
	public void addCumullatedCompound(Compound compound) {
		this.cumulatedCompounds.add(compound);
	}
	
	public void cleanCumulatedCompounds(){
		this.cumulatedCompounds.clear();
	}

	private List<Compound> byProducts = new ArrayList<Compound>();
	private List<Compound> visitedCompounds = new ArrayList<Compound>();
	private List<Compound> stopCompoundsSubstrates = new ArrayList<Compound>();
	private List<Compound> stopCompoundsByProducts = new ArrayList<Compound>();
	//private List<Compound> singleSourcesToProduceACompoundOnThePath = new ArrayList<Compound>();

	private boolean stoppedPremature = false;
	private boolean reachedSizeK = false;
	private boolean positiveNetProductionOfZeroCycle = false;
	
	boolean marked = false;
		
	private boolean stochiometricFeasibile = false;
	
	private boolean flag = false;
	private boolean emptyCompounds = false;
	private boolean containNegativeCycle = false;
	
	public boolean isMarked() {
		return marked;
	}

	public void setMarked(boolean marked) {
		this.marked = marked;
	}


	public boolean isFlag() {
		return flag;
	}

	public void setFlag(boolean flag) {
		this.flag = flag;
	}

	/**
	 * Constructors
	 *
	 */
	public PrecursorSet() {
		
		precursors = new ArrayList<Compound>();
		bootstraps = new ArrayList<Compound>();
		reactions  = new HashSet<Reaction>();
		//potentialStopCompounds = new HashSet<Compound>();
		visitedCompounds = new ArrayList<Compound>();
		setByProducts(new ArrayList<Compound>());
		//singleSourcesToProduceACompoundOnThePath = new ArrayList<Compound>();
		this.reachedSizeK = false;
		this.positiveNetProductionOfZeroCycle = false;
	}

	public PrecursorSet(int nbSources, int nbInternalCompounds, int nbReactions) {
		
		precursors = new ArrayList<Compound>();
		bootstraps = new ArrayList<Compound>();
		reactions  = new HashSet<Reaction>();
		//potentialStopCompounds = new HashSet<Compound>();
		visitedCompounds = new ArrayList<Compound>();
		setByProducts(new ArrayList<Compound>());
		//singleSourcesToProduceACompoundOnThePath = new ArrayList<Compound>();
		this.stoppedPremature = false;
		this.reachedSizeK = false;
		this.positiveNetProductionOfZeroCycle = false;

		/*this.bitSetSingleSourcesToProduceACompoundOnThePath = new byte[(int) (nbSources / 8) + 1];
		for(int i = 0; i < this.bitSetSingleSourcesToProduceACompoundOnThePath.length; i++){
			this.bitSetSingleSourcesToProduceACompoundOnThePath[i] = 0;
		}*/
	}
	
	public PrecursorSet(PrecursorSet sol) {
		
		this.precursors = new ArrayList<Compound>();
		this.bootstraps = new ArrayList<Compound>(sol.getBootstraps());
		this.reactions  = new HashSet<Reaction>();
		this.stopCompoundsByProducts = new ArrayList<Compound>();
		this.stopCompoundsSubstrates = new ArrayList<Compound>();
		//potentialStopCompounds = new HashSet<Compound>(sol.getPotentialStopCompounds());
		this.visitedCompounds = new ArrayList<Compound>();
		this.byProducts = new ArrayList<Compound>();
		//this.singleSourcesToProduceACompoundOnThePath = new ArrayList<Compound>();
		//this.bitSetSources = new byte[sol.getBitSetSources()==null? 0 : sol.getBitSetSources().length];
		for(Compound c : sol.getPrecursors()){
			addPrecursor(c);
		}
		
		// add stop comounds
		for(Compound c : sol.getStopCompoundsByProducts()){
			addByProductAsStopCompound(c);
		}
		for(Compound c : sol.getStopCompoundsSubstrates()){
			addSubstrateAsStopCompound(c);
		}
		
		// add visited compounds
		for(Compound c : sol.getVisitedCompounds()){
			addVisitedCompound(c);
		}
		
		// add by-products
		for(Compound c : sol.getByProducts()){
			addByProduct(c);
		}
		
		// add reactions
		for(Reaction r : sol.getReactions()){
			addReaction(r);
		}
		
		// if all sources can produce a compound m than we put these sources in this list
		/*this.bitSetSingleSourcesToProduceACompoundOnThePath = new byte[sol.getBitSetSingleSourcesToProduceACompoundOnThePath()==null? 0 : sol.getBitSetSingleSourcesToProduceACompoundOnThePath().length];
		for(Compound c : sol.getSingleSourcesToProduceACompoundOnThePath()){
			addSingleSourcesToProduceACompoundOnThePath(c);
		}*/
		
		this.stoppedPremature = sol.isStoppedPremature();
		this.reachedSizeK = sol.hasReachedSizeK();
		this.positiveNetProductionOfZeroCycle = sol.hasPositiveNetProductionOfZeroCycle();
	}
	
    public PrecursorSet(PrecursorSet sol, boolean stochiometricFeasibile) {
		
		this.precursors = new ArrayList<Compound>();
		this.bootstraps = new ArrayList<Compound>(sol.getBootstraps());
		this.reactions  = new HashSet<Reaction>();
		this.stopCompoundsByProducts = new ArrayList<Compound>();
		this.stopCompoundsSubstrates = new ArrayList<Compound>();
		//potentialStopCompounds = new HashSet<Compound>(sol.getPotentialStopCompounds());
		this.visitedCompounds = new ArrayList<Compound>();
		this.byProducts = new ArrayList<Compound>();
		//this.singleSourcesToProduceACompoundOnThePath = new ArrayList<Compound>();
		
		// add precursers
		for(Compound c : sol.getPrecursors()){
			addPrecursor(c);
		}
		
		// add stop comounds
		for(Compound c : sol.getStopCompoundsByProducts()){
			addByProductAsStopCompound(c);
		}
		for(Compound c : sol.getStopCompoundsSubstrates()){
			addSubstrateAsStopCompound(c);
		}
		
		// add visited compounds
		for(Compound c : sol.getVisitedCompounds()){
			addVisitedCompound(c);
		}
		
		// add by-products
		for(Compound c : sol.getByProducts()){
			addByProduct(c);
		}
		
		// add reactions
		for(Reaction r : sol.getReactions()){
			addReaction(r);
		}
		
		// if all sources can produce a compound m than we put these sources in this list
		/*this.bitSetSingleSourcesToProduceACompoundOnThePath = new byte[sol.getBitSetSingleSourcesToProduceACompoundOnThePath()==null? 0 : sol.getBitSetSingleSourcesToProduceACompoundOnThePath().length];
		for(Compound c : sol.getSingleSourcesToProduceACompoundOnThePath()){
			addSingleSourcesToProduceACompoundOnThePath(c);
		}*/
		
		this.stochiometricFeasibile = stochiometricFeasibile;
		this.stoppedPremature = sol.isStoppedPremature();
		this.reachedSizeK = sol.hasReachedSizeK();
		this.positiveNetProductionOfZeroCycle = sol.hasPositiveNetProductionOfZeroCycle();
	}

    public void union(PrecursorSet sol) {
    	for(Compound c : sol.getPrecursors()){
    		addPrecursor(c);
    	}
    	for(Compound c : sol.getBootstraps()){
    		addBootstrap(c);
    		
    	}
		
		for(Compound c : sol.getStopCompoundsByProducts()){
			addByProductAsStopCompound(c);
		}
		for(Compound c : sol.getStopCompoundsSubstrates()){
			addSubstrateAsStopCompound(c);
		}
		for(Compound c : sol.getVisitedCompounds()){
			addVisitedCompound(c);
		}
		for(Compound c : sol.getByProducts()){
			addByProduct(c);
		}
		for(Reaction r : sol.getReactions()){
			addReaction(r);
		}
		/*for(Compound c : sol.getSingleSourcesToProduceACompoundOnThePath()){
			addSingleSourcesToProduceACompoundOnThePath(c);
		}*/
		
		if(sol.isStoppedPremature()){
			this.stoppedPremature = true;
		}
		if(sol.hasReachedSizeK()){
			this.reachedSizeK = true;
		}
		if(sol.hasPositiveNetProductionOfZeroCycle()){
			this.positiveNetProductionOfZeroCycle = true;
		}
	}

	@Override
	public boolean equals(Object obj) {
		if( !(obj instanceof PrecursorSet) ) {
			return false;
		}
		PrecursorSet other = (PrecursorSet)obj;
		
		// check if the two objects have the same precursor set
		if( precursors.size() != other.precursors.size()) {
			return false;
		}
		if(haveSameSources(other) == false){
			return false;
		}
		
		// check if stop compounds are the same
		if(haveSameStopCompounds(other) == false){
			return false;
		}
				
		// check if the two objects have the same bootstrap set
		if( bootstraps.size() != other.bootstraps.size()) {
			return false;
		}
		ArrayList<Compound> myList = new ArrayList<Compound>(bootstraps);
		for(Compound c: other.bootstraps) {
			if( myList.contains(c) ) {
				myList.remove(c);
			}
		}
		if( myList.size() != 0 ) // if they are equal, I removed all element of the other list
			return false;
		
		
		//check if reactions are equal
		/*if(reactions.size() != other.reactions.size()){
			return false;
		}
		List<Reaction> myReactionList = new ArrayList<Reaction>(reactions);
		for(Reaction r : other.reactions){
			if(myReactionList.contains(r)){
				myReactionList.remove(r);
			}
		}
		if(myReactionList.size() != 0){
			return false;
		}
		*/
		return true;
	}

	
	public boolean isMinimal(PrecursorSet sol, Boolean lookBootstraps) {
		
		if(lookBootstraps == true) {
			// We minimize the union between the set of precursors and the set of bootstraps
			
			List<Compound> compoundsInA = new ArrayList<Compound>();
			compoundsInA.addAll(this.getPrecursors());
			compoundsInA.addAll(this.getBootstraps());

			List<Compound> compoundsInB = new ArrayList<Compound>();
			compoundsInB.addAll(sol.getPrecursors());
			compoundsInB.addAll(sol.getBootstraps());
			
			return compoundsInB.containsAll(compoundsInA);
			
		}
		else {
			
			if(this.getPrecursors().equals(sol.getPrecursors())) {
				// The two solutions have the same number of precursors
				// We look for the bootstrap compounds

				if(this.getBootstraps().size() == 0)
					return false;

				return sol.getBootstraps().containsAll(this.getBootstraps());
				
			}
			else {
				
				return sol.getPrecursors().containsAll(this.getPrecursors());
				
			}
		}
		
	}
	
	public Boolean isEmpty() {
		
		return(this.getPrecursors().size() == 0 && this.getBootstraps().size() == 0);
		
	}
	
	
	public void addPrecursor(Compound p) {
		if(! this.getPrecursors().contains(p)){
			this.getPrecursors().add(p);
		}
	}
	
	public void addPrecursors(List<Compound> l){
		for(Compound c : l){
			addPrecursor(c);
		}
	}
	
	public void addBootstrap(Compound b) {
		if(! this.bootstraps.contains(b)){
			this.bootstraps.add(b);
		}
	}
	
	public void add(PrecursorSet sol) {
		
		for(Compound c : sol.getPrecursors()){
    		addPrecursor(c);
    	}
    	for(Compound c : sol.getBootstraps()){
    		addBootstrap(c);
    		
    	}
		reactions.addAll(sol.getReactions());
		
		for(Compound c : sol.getStopCompoundsByProducts()){
			addByProductAsStopCompound(c);
		}
		for(Compound c : sol.getStopCompoundsSubstrates()){
			addSubstrateAsStopCompound(c);
		}
		for(Compound c : sol.getVisitedCompounds()){
			addVisitedCompound(c);
		}
		for(Compound c : sol.getByProducts()){
			addByProduct(c);
		}
	}

	public List<Compound> getBootstraps() {
		return bootstraps;
	}

	public void setBootstraps(List<Compound> bootstraps) {
		this.bootstraps = bootstraps;
	}

	public List<Compound> getPrecursors() {
		return precursors;
	}

	public void setPrecursors(List<Compound> precursors) {
		removeAllPrecursors();
		for(Compound p : precursors){
			addPrecursor(p);
		}
	}
	
	public Set<Reaction> getReactions() {
		return reactions;
	}

	public void setReactions(Set<Reaction> reactions) {
		this.reactions = reactions;
	}

	public String toString() {
		
		/*HashSet <String> compounds = new HashSet<String>();
		for(Reaction r : reactions){
			compounds.addAll(r.getSubstrates().keySet());
			compounds.addAll(r.getProduces().keySet());
		}*/
		
		
		//String xml = "Precursors: "+precursors + " Reactions: " + reactions + " Stop compounds: [" + getStopCompoundsByProducts() + ", " + getStopCompoundsSubstrates() + "]";
		String xml = "Precursors: "+precursors + " Stop compounds: [" + getStopCompoundsByProducts() + ", " + getStopCompoundsSubstrates() + "]";
		
		
		/*
		xml = xml + "\n\n<sbml xmlns=\"http://www.sbml.org/sbml/level2\" version=\"1\" level=\"2\" xmlns:html=\"http://www.w3.org/1999/xhtml\"><model id=\"test\">";
		xml = xml + "<listOfSpecies>";
		for(String c : compounds){
			xml = xml + "<species id=\"" + c + "\" name=\"" + c + "\"/>";
		}
		xml = xml + "</listOfSpecies><listOfReactions>";
		for(Reaction r : reactions){
			xml = xml + "<reaction id=\"" + r.getId() +"\" name=\"" + r.getName() + "\" reversible=\"" + r.isReversible() + "\">";
			Iterator<String> iter = r.getSubstratesStoich().keySet().iterator();
			xml = xml + "<listOfReactants>";
			while(iter.hasNext()){
				String c = iter.next();
				xml = xml + "<speciesReference species=\"" + c + "\" stoichiometry=\"" + r.getSubstratesStoich().get(c) + "\"/>";
			}
			xml = xml + "</listOfReactants>";
			xml = xml + "<listOfProducts>";
			iter = r.getProductsStoich().keySet().iterator();
			while(iter.hasNext()){
				String c = iter.next();
				xml = xml + "<speciesReference species=\"" + c + "\" stoichiometry=\"" + r.getProductsStoich().get(c) + "\"/>";
			}
			xml = xml + "</listOfProducts>";
			xml = xml + "</reaction>";
		}
		xml = xml + "</listOfReactions></model></sbml>";
		*/
		return xml;
//		return "Precursors: "+precursors + " Reactions: " + reactions;
	}

	@Override
	public int compareTo(PrecursorSet ps) {
		return toString().compareTo(ps.toString());
		//return this.hashCode() - ps.hashCode();
	}

	public boolean isEmptyCompounds() {
		return emptyCompounds;
	}

	public void setEmptyCompounds(boolean emptyCompounds) {
		this.emptyCompounds = emptyCompounds;
	}

	/*
	 * Sets in arbitrary form the stochiometricFeasible flag.
	 * Automatically sets this PrecursorSet as changed, in order to allow
	 * that next call to checkStoichiometricFeasibility works properly.
	 */
	public void setStochiometricFeasibile(boolean stochiometricFeasibile) {
		this.stochiometricFeasibile = stochiometricFeasibile;
	}
	
	public boolean isStochiometricFeasibile() {
		return stochiometricFeasibile;
	}

	public void addReaction(Reaction r) {
		if(! this.reactions.contains(r)){
			this.reactions.add(r);
			//this.changedAfterFeasibilityCheck = true;
		}	
	}
	
	public void addReactions(Set<Reaction> s) {
		for(Reaction r : s){
			addReaction(r);
		}
	}
	
	public Set<Compound> getAllImpliedCompounds(Compound target){
		HashSet<Compound> returnSet = new HashSet<Compound>();
		int countTargetAsProduct = 0;
		for (Reaction r: reactions){
			returnSet.addAll(r.getSubstrates().values());
			for(Compound c : r.getProduces().values()){
				returnSet.add(c);
				if(target != null && c.getId().equals(target.getId())){
					countTargetAsProduct++;
				}
			}
		}
		if(countTargetAsProduct == 1){
			Iterator<Compound> iter = returnSet.iterator();
			while(iter.hasNext()){
				Compound c = iter.next();
				if(target != null && c.getId().equals(target.getId())){
					iter.remove();
				}
			}
		}
		return returnSet;
	}
	
	public Set<Compound> getAllSideProducts(Compound target){
		HashSet<Compound> byProducts = new HashSet<Compound>();
/*		int countTargetAsProduct = 0;
		for(Reaction r : reactions){
			//byProducts.addAll(r.getProduces().values());
			for(Compound c : r.getProduces().values()){
				if(c.getId().equals("TARGET_MINIMAL") || c.getId().equals("TARGET")){
					byProducts.remove(c);
				}
				else if(c.getId().equals(target.getId())){
					countTargetAsProduct++;
					byProducts.add(c);
				}
				else{
					byProducts.add(c);
				}
			}
		}
		for(Reaction r : reactions){
			for(Compound c : r.getSubstrates().values()){
				if(byProducts.contains(c)){
					byProducts.remove(c);
				}
			}
		}
		
		// if the target in scope appears more than one time in the reaction set
		// than it is also a by-product
		if(countTargetAsProduct == 1){
			byProducts.remove(target);
		}
		
		*/
	
		HashMap<Compound, Double>  newAmount = new HashMap<Compound, Double>();
		Iterator<Reaction> iterReac = this.reactions.iterator();
		double fluxFactor = 1.0;
		int count = 0;
		Reaction lastReaction = null;
		while(iterReac.hasNext()){
			count++;
			Reaction r = iterReac.next();
			double newFluxFactor = 1.0;
			if(count > 1){
				//target = ...;
				newFluxFactor = fluxFactor * lastReaction.getSubstratesStoich().get(target.getId()) / r.getProductsStoich().get(target.getId());
			}
			for(Compound c: r.getProduces().values()){
				//if(! c.getId().equals(a.getId())){ // just for side products -> compute new amount
					if(newAmount.containsKey(c.getId())){
						//System.out.println("AM P " + c.getId() + ": " + newAmount.get(c.getId()) + " + " + r.getProductsStoich().get(c.getId()) + " * " + factor);
						newAmount.put(c, newAmount.get(c.getId()) + (newFluxFactor * r.getProductsStoich().get(c.getId())) );
					}
					else if(! c.getId().equals(target.getId())){
						//System.out.println("AM P " + c.getId() + ": "  + r.getProductsStoich().get(c.getId()) + " * " + factor);
						newAmount.put(c, (newFluxFactor * r.getProductsStoich().get(c.getId())) );
					}
			//}
			}
		
			for(Compound c: r.getSubstrates().values()){
				if(newAmount.containsKey(c.getId())){
					//System.out.println("AM S " + c.getId() + ": " + newAmount.get(c.getId()) + " - " + r.getSubstratesStoich().get(c.getId()) + " * " + factor);
					newAmount.put(c, newAmount.get(c.getId()) - (newFluxFactor * r.getSubstratesStoich().get(c.getId())) );
				}
			}
			
			fluxFactor = newFluxFactor;
			lastReaction = r;
		}
		
		
		return byProducts;
	}
	
	static class PropertiesReturnType{
		boolean allLeafesAreSources;
		boolean hasCycle;
		
		PropertiesReturnType(boolean allLeafesAreSources, boolean hasCycle){
			this.allLeafesAreSources = allLeafesAreSources;
			this.hasCycle = hasCycle;
		}
	}

	/*
	 * We never call this function
	 */
	/*public boolean combinableWith(PrecursorSet c2, Compound target) {
	 
		// TODO Auto-generated method stub
		
		Set<Compound> c1Compounds = this.getAllImpliedCompounds(target);
		Set<Compound> c2Compounds = c2.getAllImpliedCompounds(target);
		
		boolean sideProductInPS2 = false;
		boolean sideProductInPS1 = false;
		if(this.getAmountOfCompounds().size() == 0){
			return false;
		}
		if(c2.getAmountOfCompounds().size() == 0){
			return false;
		}
		for(String cId : this.getAmountOfCompounds().keySet()){
			if(this.getAmountOfCompounds().get(cId) <= 0.0){
				continue;
			}
			for(Compound c : c2Compounds){
				if(cId.equals(c.getId())){
					sideProductInPS2 = true;
					break;
				}
			}
			if(sideProductInPS2){
				break;
			}
		}
		
		
		for(String cId : c2.getAmountOfCompounds().keySet()){
			if(c2.getAmountOfCompounds().get(cId) <= 0.0){
				continue;
			}
			for(Compound c : c1Compounds){
				if(cId.equals(c.getId())){
					sideProductInPS1 = true;
					break;
				}
			}
			if(sideProductInPS1){
				break;
			}
		}
		
		
		if(sideProductInPS1 && sideProductInPS2){
			return true;
		}
		return false;
	}
*/
	public List<Compound> getStopCompounds() {
		List<Compound> allStopCompounds = new ArrayList<Compound>();
		allStopCompounds.addAll(this.stopCompoundsByProducts);
		allStopCompounds.addAll(this.stopCompoundsSubstrates);
		return allStopCompounds;
	}

	public void setByProductsAsStopCompounds(List<Compound> stopCompounds) {
		removeAllByProductsAsStopCompounds();
		for(Compound s : stopCompounds){
			addByProductAsStopCompound(s);
		}
	}
	
	public void setSubstratesAsStopCompounds(List<Compound> stopCompounds) {
		removeAllSubstratesAsStopCompounds();
		for(Compound s : stopCompounds){
			addSubstrateAsStopCompound(s);
		}
	}
	
	public void addByProductAsStopCompound(Compound c){
		if(! this.stopCompoundsByProducts.contains(c)){
			this.stopCompoundsByProducts.add(c);
		}
	}
	
	public void addSubstrateAsStopCompound(Compound c){
		if(! this.stopCompoundsSubstrates.contains(c)){
			this.stopCompoundsSubstrates.add(c);
			}
	}
	
	/*public void removeStopCompound(Compound c){
		Iterator<Compound> iter = this.stopCompounds.iterator();
		while(iter.hasNext()){
			if(iter.next().equals(c)){
				iter.remove();
				if(this.bitSetStopCompounds != null){
					int idxByte = (int) c.getIdxBitSet() / 8;
					int idxBit = (int) c.getIdxBitSet() % 8;
					this.bitSetStopCompounds[idxByte] &= ~(1 << idxBit);
				}
				this.changedAfterFeasibilityCheck = true;
			}
		}
		
	}*/
	
	public void removeByProductAsStopCompound(Compound c){
		Iterator<Compound> iter = this.stopCompoundsByProducts.iterator();
		while(iter.hasNext()){
			if(iter.next().equals(c)){
				iter.remove();
				break;
			}
		}
	}

	public void removeSubstrateAsStopCompound(Compound c){
		Iterator<Compound> iter = this.stopCompoundsSubstrates.iterator();
		while(iter.hasNext()){
			if(iter.next().equals(c)){
				iter.remove();
				break;
			}
		}
	}
	
	public void removeVisitedCompound(Compound c){
		Iterator<Compound> iter = this.visitedCompounds.iterator();
		while(iter.hasNext()){
			if(iter.next().equals(c)){
				iter.remove();
				break;
			}
		}
	}

	public void removeAllByProductsAsStopCompounds(){
		this.stopCompoundsByProducts.clear();
	}
	
	
	public void removeAllSubstratesAsStopCompounds(){
		this.stopCompoundsSubstrates.clear();		
	}
	
	public void removeAllPrecursors(){
		for(Compound c : this.precursors){
			removePrecursor(c);
		}
		this.precursors.clear();
	}
	
	public void removePrecursor(Compound c){
		this.precursors.remove(c);
	}

	public boolean containNegativeCycle() {
		return containNegativeCycle;
	}

	public void setNegativeCycle(boolean containNegativeCycle) {
		this.containNegativeCycle = containNegativeCycle;
	}

/*	public Set<Compound> getPotentialStopCompounds() {
		return potentialStopCompounds;
	}

	public void setPotentialStopCompounds(Set<Compound> potentialStopCompounds) {
		this.potentialStopCompounds = potentialStopCompounds;
	}
	
	public void addPotentialStopCompounds(Compound c){
		this.potentialStopCompounds.add(c);
	}
*/
	
	public List<Compound> getVisitedCompounds() {
		return visitedCompounds;
	}

	public void setVisitedCompounds(List<Compound> visitedCompounds) {
		removeAllVisitedCompounds();
		for(Compound s : visitedCompounds){
			addVisitedCompound(s);
		}
	}
	
	private void removeAllVisitedCompounds() {
		// TODO Auto-generated method stub
		this.visitedCompounds.clear();
	}

	public void addVisitedCompound(Compound c){
		if(! this.visitedCompounds.contains(c)){
			this.visitedCompounds.add(c);
		}
	}
	
	public void addVisitedCompounds(Set<Compound> s){
		for(Compound c : s){
			addVisitedCompound(c);
		}
	}

	public List<Compound> getByProducts() {
		return byProducts;
	}

	public void setByProducts(List<Compound> byProducts) {
		removeAllByProducts();
		for(Compound s : byProducts){
			addByProduct(s);
		}
	}
	
	private void removeAllByProducts() {
		// TODO Auto-generated method stub
		this.byProducts.clear();
	}

	public void addByProduct(Compound c){
		if(! this.byProducts.contains(c) && ! c.isPrecursor()){
			this.byProducts.add(c);
		}
	}

	/*public void addSingleSourcesToProduceACompoundOnThePath(Compound c){
		if(! this.singleSourcesToProduceACompoundOnThePath.contains(c)){
			this.singleSourcesToProduceACompoundOnThePath.add(c);
			if(this.bitSetSingleSourcesToProduceACompoundOnThePath != null){
				int idxByte = (int) c.getIdxBitSet() / 8;
				int idxBit = (int) c.getIdxBitSet() % 8;
				this.bitSetSingleSourcesToProduceACompoundOnThePath[idxByte] |= (1 << idxBit);
			}
		}
	}*/
	
	public void addByProducts(Collection<Compound> l){
		for(Compound c : l){
			addByProduct(c);
		}
	}
	
	public void removeByProduct(Compound c){
		this.byProducts.remove(c);
	}

	public boolean haveSameSources(PrecursorSet ps){
		if(this.getPrecursors().size() != ps.getPrecursors().size()){
			return false;
		}
		HashSet<Compound> hashSources = new HashSet<Compound>();
		for(int i = 0; i < ps.getPrecursors().size(); i++)
			hashSources.add(ps.getPrecursors().get(i));
		for(int j = 0; j < this.getPrecursors().size(); j++)
		{
			if( ! hashSources.contains(this.getPrecursors().get(j)) )
				return false;
		}
		return true;
	}
	
	public boolean haveSameStopCompounds(PrecursorSet ps){
		HashSet<Compound> hashStopCompounds = new HashSet<Compound>();
		for(int i = 0; i < ps.getStopCompounds().size(); i++)
			hashStopCompounds.add(ps.getStopCompounds().get(i));
		for(int j = 0; j < this.getStopCompounds().size(); j++)
		{
			if( ! hashStopCompounds.contains(this.getStopCompounds().get(j)) )
				return false;
		}
		
		return (this.getStopCompounds().size() == ps.getStopCompounds().size());
	}
	
	
	public boolean isStopCompoundSubSetOf(PrecursorSet ps, boolean properSubset){
		HashSet<Compound> hashStopCompounds = new HashSet<Compound>();
		for(int i = 0; i < ps.getStopCompounds().size(); i++)
			hashStopCompounds.add(ps.getStopCompounds().get(i));
		for(int j = 0; j < this.getStopCompounds().size(); j++)
		{
			if( ! hashStopCompounds.contains(this.getStopCompounds().get(j)) )
				return false;
		}
		
		return !properSubset || (this.getStopCompounds().size() < ps.getStopCompounds().size());
	}
	
	public boolean isSourcesSubSetOf(PrecursorSet ps, boolean properSubset){
		HashSet<Compound> hashSources = new HashSet<Compound>();
		for(int i = 0; i < ps.getPrecursors().size(); i++)
			hashSources.add(ps.getPrecursors().get(i));
		for(int j = 0; j < this.getPrecursors().size(); j++)
		{
			if( ! hashSources.contains(this.getPrecursors().get(j)) )
				return false;
		}
		
		return !properSubset || (this.getPrecursors().size() < ps.getPrecursors().size());
	}
	
	public boolean isReactionSubSetOf(PrecursorSet ps, boolean properSubset){
		for(Reaction r : this.getReactions())
		{
			if( ! ps.getReactions().contains(r) )
				return false;
		}
		
		return !properSubset || (this.getReactions().size() < ps.getReactions().size());
	}
	
	public Set<Reaction> unionOfReactions(PrecursorSet ps){
		Set<Reaction> unionOfReactions = this.getReactions();
		unionOfReactions.addAll(ps.getReactions());
		return unionOfReactions;
	}

	public boolean hasIntersectionBetweenByProductAndLineCompound(
			PrecursorSet ps) {
		// TODO Auto-generated method stub
		boolean intersection = false;
		for(Compound byProduct : this.byProducts){
			for(Compound visitedCompound : ps.getVisitedCompounds()){
				if(byProduct.equals(visitedCompound)){
					intersection = true;
					break;
				}
			}
			if(intersection == true){
				break;
			}
		}
		if(intersection == false){
			return false;
		}
		// opposite direction
		for(Compound byProduct : ps.byProducts){
			for(Compound visitedCompound : this.getVisitedCompounds()){
				if(byProduct.equals(visitedCompound)){
					return intersection;
				}
			}
		}
		return false;
	}

	public List<Compound> getStopCompoundsSubstrates() {
		return stopCompoundsSubstrates;
	}

	public List<Compound> getStopCompoundsByProducts() {
		return stopCompoundsByProducts;
	}

	public void addSubstrateAsStopCompounds(
			List<Compound> notYetVisitedNonSourceCompounds) {
		// TODO Auto-generated method stub
		for(Compound c : notYetVisitedNonSourceCompounds){
			addSubstrateAsStopCompound(c);
		}
	}

	public boolean isStoppedPremature() {
		return stoppedPremature;
	}

	public void setStoppedPremature(boolean stoppedPremature) {
		this.stoppedPremature = stoppedPremature;
	}

	public boolean hasReachedSizeK() {
		return reachedSizeK;
	}

	public void setReachedSizeK(boolean reachedSizeK) {
		this.reachedSizeK = reachedSizeK;
	}

	public boolean hasPositiveNetProductionOfZeroCycle() {
		return positiveNetProductionOfZeroCycle;
	}

	public void setPositiveNetProductionOfZeroCycle(
			boolean positiveNetProductionOfZeroCycle) {
		this.positiveNetProductionOfZeroCycle = positiveNetProductionOfZeroCycle;
	}

	/*public List<Compound> getSingleSourcesToProduceACompoundOnThePath() {
		return singleSourcesToProduceACompoundOnThePath;
	}

	public void setSingleSourcesToProduceACompoundOnThePath(
			List<Compound> singleSourcesToProduceACompoundOnThePath) {
		this.singleSourcesToProduceACompoundOnThePath = singleSourcesToProduceACompoundOnThePath;
	}

	public void addSingleSourcesToProduceACompoundOnThePath(
			List<Compound> sources) {
		for(Compound c : sources){
			addSingleSourcesToProduceACompoundOnThePath(c);
		}
	}
*/
	
}
