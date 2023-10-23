/*SASITA - A software tool to find all minimal precursor sets for a given set of targets in metabolic networks.
Copyleft (C) 2015 Ricardo Andrade (ricardoluizandrade@gmail.com), Martin Wannagat (mwannagat@gmail.com)
Copyleft (C) 2011 Ludovic Cottret (l.cottret@gmail.com), Paulo Vieira Milreu  (paulovieira@milreu.com.br)      
This file is part of SASITA.

SASITA is free software: you can redistribute it and/or modify
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
import java.util.Collections;
import java.util.Date;
import java.util.List;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import metabolicNetwork.Compound;
import metabolicNetwork.MetabolicNetwork;
import metabolicNetwork.Reaction;
import pitufolandia.Pitufrankstein;
import pitufolandia.Sasita;
import utils.ArgumentParser;
import utils.SolutionChecker;
import utils.StringUtils;

public class Main {

	MetabolicNetwork network;
	MetabolicNetwork fullNetwork;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) 
	{
		ArgumentParser p = new ArgumentParser(args);
		if (args.length == 0 || p.hasOption("h")) 
			printUsage();
		else 
		{
			if(! p.hasOption("s") || !p.hasOption("i")) 
				printUsage();
			System.out.println("Started at "+new Date());
			
			InputParameters.oneByOne = p.hasOption("o");
			InputParameters.nbThreads = p.hasOption("nbThreads") ? Integer.parseInt(p.getOption("nbThreads")) : 1;
			InputParameters.considerOnlyUserDefinedPrecursors =  p.hasOption("udpo") ? true : false;
			InputParameters.noPreprocessing = p.hasOption("nopreprocess") ? true : false;
			InputParameters.verbose =  p.hasOption("verbose") ? true : false;
			InputParameters.checkMinimality = p.hasOption("checkMin");
			InputParameters.fbaCheck = p.hasOption("fbaCheck") ? true : false;
			InputParameters.epsilon1 = p.hasOption("epsilon1") ? Double.parseDouble(p.getOption("epsilon1")) : 0.1;
			InputParameters.epsilon2 = p.hasOption("epsilon2") ? Double.parseDouble(p.getOption("epsilon2")) : 0.1;
			InputParameters.bigM = p.hasOption("bigM") ? Double.parseDouble(p.getOption("bigM")) : 1000;
			InputParameters.modeDuplicatingMachinery = p.hasOption("dupMach") ? true : false;
			InputParameters.modeSteadyState = p.hasOption("steadyState") ? true : false;
			InputParameters.modeCheckSolutionsOnly = p.hasOption("checkSolsOnly") ? true : false;
			InputParameters.addReactionsToSolutions = p.hasOption("addReactions") ? true : false;
			InputParameters.addCumulatedCompoundsToSolutions = p.hasOption("addCumulated") ? true : false;
			InputParameters.TiLim = p.hasOption("tiLim") ? Integer.parseInt(p.getOption("tiLim")) : 3600;


			String sbmlFile = p.getOption("s");
			String inputFile = p.getOption("i");
			String solutionsToCheckFile = p.getOption("p");
			
			InputParameters.filterPairedCofactors = p.hasOption("fCof") ? p.getOption("fCof") : "";

			if (InputParameters.fbaCheck && (InputParameters.modeDuplicatingMachinery || InputParameters.modeSteadyState)){
				System.err.println("[ERROR] -fbaCheck does not work with -dupMach or -steadyState.");
				printUsage();
			}
			
			if (InputParameters.modeDuplicatingMachinery && InputParameters.modeSteadyState){
				System.err.println("[ERROR] Incompatible modes. You should chose either -dupMach or -steadyState.");
				printUsage();
			}
			
			System.out.println("Number of threads: " + InputParameters.nbThreads);

			printGPLShortNotice();
			System.out.println("SBML File:"+sbmlFile);
			System.out.println("Input File:"+inputFile);

			Main main = new Main(sbmlFile, inputFile);

			// Extend topological precursors (mark as precursor), bootstraps, do the forward propagation step, the maximal target reduction
			
			main.preprocessNetwork();
			
			// Check if in mode of checking solutions only
			if (InputParameters.modeCheckSolutionsOnly){
				try {
					SolutionChecker checker = new SolutionChecker(main.getNetwork(), InputParameters.getTargetCompounds());
					checker.checkSolutions(solutionsToCheckFile);
				} catch (Exception e) {
					System.err.println(e.getMessage());
				}
				System.out.println("\nFinished at "+new Date());
				System.exit(0);
			}
			
			// Launch SASITA
			
			long startTime = System.nanoTime();    
			
			Sasita finder = new Sasita(main.getNetwork(), false);
			finder.setEpsilon1(InputParameters.epsilon1);
			finder.setEpsilon2(InputParameters.epsilon2);
			finder.setBigM(InputParameters.bigM);
			finder.setDuplicatingMachineryActive(InputParameters.modeDuplicatingMachinery);
			finder.setSteadyState(InputParameters.modeSteadyState);
			List<PrecursorSet> solutions = finder.findPrecursorsInNetwork(InputParameters.getTargetCompounds());
			
			 /*Pitufrankstein finderTest = new Pitufrankstein(main.getNetwork(), false);
			 finderTest.setEpsilon1(InputParameters.epsilon1);
			 finderTest.setBigM(InputParameters.bigM);
			 List<PrecursorSet> solutions = finderTest.findPrecursorsInNetwork(InputParameters.getTargetCompounds());
			*/
			if (solutions != null){
				finder.printSolutions2Xml(solutions, sbmlFile + ".output");
			}
			if (InputParameters.checkMinimality){
				finder.checkMinimality(solutions);
			}
			
			long estimatedTime = System.nanoTime() - startTime;
			
			System.out.println("\nElapsed time: " + TimeUnit.NANOSECONDS.toMillis(estimatedTime) + "ms");
			System.out.println("\nFinished at "+new Date());
		}
	}
	public static void printGPLShortNotice()
	{
		System.out.println("------------------------------SASITA------------------------------");
		System.out.println("Method coded by Andrade & Wannagat");
		System.out.println("Some modules coded by Cottret & Milreu");
		System.out.println("This is a free software, under CeCILL license");
		System.out.println("you are welcome toredistribute it,");
		System.out.println("under certain conditions.");
		System.out.println("Check the COPYING file on the sources directory ");
		System.out.println("for more details.");
	}	

	public static void printUsage() {
		printGPLShortNotice();
		System.err.println("-----------------------SASITA Main options-----------------------");
		System.err.println("-s=filename\tSBML file (mandatory)");
		System.err.println("-i=filename\tInput file (mandatory)");
		System.err.println("-o\t\tLook separately at each target compound");
		System.err.println("-epsilon1=X\tSet the value of epsilon1 as X (default: 0.1)");
		System.err.println("-epsilon2=X\tSet the value of epsilon2 as X (default: 0.1)");
		System.err.println("-bigM=X\t\tSet the value of bigM as X (default: 1000)");
		System.err.println("-udpo\t\tConsider only user defined precursors as source\n\t\tcompounds");
		System.err.println("-nbThreads=X\tSet as X the amount of threads for parallel \n\t\tprocessing (default: 1)");
		System.err.println("-nopreprocess\tSkip optional preprocessing");
		System.err.println("-addReactions\tAdd the IDs of the trigered reactions to each \n\t\tsolution");
		System.err.println("-addCumulated\tAdd the IDs of the compounds that accumulate (Sv>0)\n\t to each \n\t\tsolution");
		System.err.println("-checkMin\tCheck if all solutions are minimal (helps \n\t\tdetecting if parameters were correctly choosed)");
		System.err.println("-fbaCheck\tCheck if each solution can really produce the\n\t\ttarget through a FBA test where all compounds can \n\t\taccumulate (helps detecting if parameters were \n\t\tcorrectly choosed)");
		System.err.println("-verbose\tPrints more information in the screen");
		System.err.println("--------------------------Optional Modes--------------------------");
		System.err.println("You can choose (optionally) one of the following modes:");
		System.err.println("-dupMach\tRun using the duplicating machinery model");
		System.err.println("-steadyState\tRun using the steady state model (Sv=0)\n");
		System.err.println("-checkSolsOnly\tRun just a check of solutions given in the XML \n\t\tformat through the flag \"-p=<file>\" \n");
		System.err.println("If no optional mode is selected, the default mode is ");
		System.err.println("the Stoichiometric Factory mode (Sv>=0).");
		System.err.println("------------------------------------------------------------------");
		System.exit(-1);
	}

	public Main(MetabolicNetwork network, Set<String> inputIds, Set<String> bootstrapIds, Set<String> targetIds, Set<String> precursorIds) {

		this.setNetwork(network);

		for(String id : inputIds) {
			Compound c = this.getNetwork().getCompounds().get( id );
			if( c == null)
				System.out.println("Input compound " + id +" not in the network.");
			else
			{
				InputParameters.getInputCompounds().add(c);
			}
		}

		for(String id : bootstrapIds) {
			Compound c = this.getNetwork().getCompounds().get( id );
			if( c == null)
				System.out.println("Bootstrap Compound " + id + " not in the network.");
			else
			{
				InputParameters.getBootstrapCompounds().add(c);
			}
		}

		for(String id : targetIds) {
			Compound c = this.getNetwork().getCompounds().get( id );
			if( c == null) {
				System.out.println("Target Compound " + id + " not in the network.");
				// Creation
				c = new Compound(id, id, "");
				this.getNetwork().getCompounds().put(id, c);
			}
			else
			{
				InputParameters.getTargetCompounds().add(c);
			}
		}

		for(String id : precursorIds) {
			Compound c = this.getNetwork().getCompounds().get( id );
			if( c == null)
				System.out.println("Precursor Compound " + id + " not in the network.");
			else
				c.setUserDefinedPrecursor(true);
		}
	}


	public Main(String sbmlFile, String inputFile) 
	{
			MetabolicNetwork network = new MetabolicNetwork();		
			network.parseSbmlFormat(sbmlFile);
			
			System.out.println("\nThe original network contains "+ network.getReactions().size()+" reactions and "+network.getCompounds().size()+" compounds.");
			
			this.setNetwork(network);
			
			this.readInputFile(inputFile);
	}

	public MetabolicNetwork getNetwork() {
		return network;
	}
	
	public MetabolicNetwork getFullNetwork() {
		return fullNetwork;
	}
	

	public void preprocessNetwork()
	{	
		// mark compounds defined as bootstraps
		for(Compound cpd : InputParameters.getBootstrapCompounds()) {
			cpd.setBootstrap(true);
		}

		// remove catalyst compounds (present as substrate and product of the same reaction),
		// orphans (compounds that are used by no reaction) and "empty" reaction (reaction with empty substrate or product list)
		network.removeCompoundsFromReactionsIfPresentAsSubstrateAndProduct();
		if (!InputParameters.modeSteadyState){
			// We don't remove empty reactions in the case of the steady state.
			// These reactions are important in this case.
			network.cleanEmptyReactions();
		}
		else{
			//In this case we will remove only the reactions without substrate
			network.cleanReactionsWithNoSubstrate();
		}
		network.cleanRepeatedReactions();
		network.cleanOrphanCompounds();
		
		// remove unbalanced reactions
		if(! InputParameters.noPreprocessing){
			network.removeUnBalancedReactions();
			// remove co-factors if reaction remains still balanced
			//fullNetwork = network.backupHardCopy();
			network.removePairedCoFactors(true);
			network.removeCO2();
			network.removeSubReactions();
		}
		// first decide which are the topological precursors
		network.extendTopologicalPrecursors( InputParameters.precursorIfProducedOnlyByReversible );				

		// All compounds defined as precursors have to be transformed in topological precursors
		network.extendNonTopologicalPrecursors();

		List<Compound> listPrecursorOld = new ArrayList<Compound>();
		List<Compound> listOfbootstraps = InputParameters.getBootstrapCompounds();
		List<Compound> compoundList = new ArrayList<Compound>();
		for(Compound c : this.network.getCompounds().values()){
			if(c.isPrecursor()){
				listPrecursorOld.add(c);
			}
			compoundList.add(c);
		}
		
		// for Azrael: check if we can produce all targets if we use all sources and if we can produce something out of nothing
		if(! InputParameters.noPreprocessing){
			network.checkNetwork(InputParameters.getTargetCompounds(), compoundList);
			network.cleanOrphanCompounds();
			network.extendTopologicalPrecursors( InputParameters.precursorIfProducedOnlyByReversible );
		}

		List<Compound> listPrecursorNew = new ArrayList<Compound>();
		for(Compound c : this.network.getCompounds().values()){
			if(c.isPrecursor()){
				listPrecursorNew.add(c);
			}
		}
		// Remove reactions outside paths from sources (and bootstraps) to targets
		if(! InputParameters.noPreprocessing){
			List<Compound> sourcesAndBootstraps = new ArrayList<Compound>(listPrecursorOld);
			sourcesAndBootstraps.addAll(listOfbootstraps);
			deleteRecursivelyReactionWithUnAcceptedSources(sourcesAndBootstraps);
		}
		int nRev=0;
		for(Reaction rxn : network.getReactions().values()) 
		{
			if(rxn.isReversible())
				nRev++;
		}
		System.out.println("\nNumber of reversible reactions : "+(nRev/2)+"\n");
	}
	
	public void assignIdx2CompoundsAndReactions(){
		network.extendTopologicalPrecursors( InputParameters.precursorIfProducedOnlyByReversible );
		// asign index to each compound (precursor, other source)
		int idxSource = 0;
		int idxIntCompound = 0;
		for(Compound c : network.getCompounds().values()){
			if(c.isPrecursor()){
				c.setIdxBitSet(idxSource++);
			}
			else{
				c.setIdxBitSet(idxIntCompound++);
			}
		}
		
		// asign index to each reaction
		int idxReaction = 0;
		for(Reaction r : network.getReactions().values()){
			r.setIdxBitSet(idxReaction++);
		}
		
		
		
		
		
		if(InputParameters.verbose){
			for(Reaction r : network.getReactions().values()){
				System.out.println("Remains: " + r);
			}
		}
	}

	private void readInputFile(String inputFile)	
	{
		InputParametersParser.parse(inputFile);

		// adding the input compounds for the forward propagation process
		for(String strId : InputParametersParser.input ) {
			String id = StringUtils.sbmlEncode(strId);
			Compound c = this.getNetwork().getCompounds().get( id );
			if( c == null)
				System.out.println("Input Compound " + id + " not in the network.");
			else
			{
				InputParameters.getInputCompounds().add(c);
			}
			System.out.println("Input compound: " + strId + " - SBML Format: " + id);         
		}

		// adding the bootstrap compounds for the forward propagation process
		for(String strId : InputParametersParser.bootstrap ) {
			String id = StringUtils.sbmlEncode(strId);
			Compound c = this.getNetwork().getCompounds().get( id );
			if( c == null)
				System.out.println("Bootstrap Compound " + id + " not in the network.");
			else
			{
				InputParameters.getBootstrapCompounds().add(c);
			}
			System.out.println("Bootstrap compound: " + strId + " - SBML Format: " + id);         
		}

		// adding the user defined precursors 
		for(String strId : InputParametersParser.userDefinedPrecursor ) {
			String id = StringUtils.sbmlEncode(strId);
			Compound c = this.getNetwork().getCompounds().get( id );
			if( c == null)
				System.out.println("UD Precursor " + id + " not in the network.");
			else
			{
				c.setUserDefinedPrecursor(true);
				System.out.println("UD Precursor: " + strId + " - SBML Format: " + id);
			}
			         
		}		

		// adding the forbidden precursors
		for(String strId : InputParametersParser.forbiddenPrecursors ) {
			String id = StringUtils.sbmlEncode(strId);
			Compound c = this.getNetwork().getCompounds().get( id );
			if( c == null)
				System.out.println("Forbidden Precursor " + id + " not in the network.");
			else
			{
				c.setAllowed(false);
			}
			System.out.println("Forbidden Precursor: " + strId + " - SBML Format: " + id);         
		}		

		// adding the target precursors
		for(String strId : InputParametersParser.target ) {
			String id = StringUtils.sbmlEncode(strId);
			Compound c = this.getNetwork().getCompounds().get( id );
			if( c == null)
				System.out.println("Target " + id + " not in the network.");
			else
			{
				InputParameters.getTargetCompounds().add(c);
				System.out.println("Target: " + strId + " - SBML Format: " + id);   
			}
			      
		}		
	}

	// Removes from the network all the compounds in the list of "not allowed precursors"
	public void removeForbiddenPrecursors()
	{
		List<Compound> cpdsToRemove = new ArrayList<Compound>();
		for(Compound c : network.getCompounds().values() )
		{
			if( c.isPrecursor() && !c.isAllowed() )
				cpdsToRemove.add(c);
		}
		network.removeCompounds(cpdsToRemove);
	}

	public void setNetwork(MetabolicNetwork network) {
		this.network = network;
	}
	
	public void setFullNetwork(MetabolicNetwork network) {
		this.fullNetwork = network;
	}
	

	/**
	 * From a set of ids, returns the list of corresponding compounds in the network
	 * @param ids
	 * @return
	 */
	public List<Compound> strToCompound(String ids[], boolean sorting) {
		List<Compound> cpds = new ArrayList<Compound>();

		for(String id : ids) {
			Compound c = this.getNetwork().getCompounds().get( id );
			if( c == null)
				System.out.println("Compound " + id + " not in the network.");
			else
			{
				cpds.add(c);
			}
		}

		if(sorting) {
			Collections.sort(cpds, null);
		}

		return cpds;
	}
	
	private void deleteRecursivelyReactionWithUnAcceptedSources(List<Compound> acceptedSources){
		boolean removedReaction = false;
		do{
			removedReaction = false;
			for(Compound cmp : this.network.getCompounds().values()){
				//if it is not an allowed source
				if(!acceptedSources.contains(cmp) &&
						(cmp.getProducedBy().size() == 0 || 
						(cmp.getProducedBy().size() == 1 && cmp.getSubstrateOf().size() == 1 && cmp.getProducedBy().get(0).isReversible()))){
						List<Reaction> consumedBy = new ArrayList<Reaction>();
						consumedBy.addAll(cmp.getSubstrateOf());
						for(Reaction r : consumedBy){
							System.out.println("Remove reaction that consume " + cmp + ": " + r);
							this.network.removeReaction(r);
							removedReaction = true;
						}
				}
			}
		} while(removedReaction);
	}
}
