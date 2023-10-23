package pitufolandia;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;

import metabolicNetwork.Compound;
import metabolicNetwork.MetabolicNetwork;
import metabolicNetwork.Reaction;
import pitufo.PrecursorFinder;
import utils.CPLEXInterface;
import utils.MetabolicNetworkSBMLWriter;
import utils.OptimisationInterface;
import utils.SolutionChecker;
import utils.StringUtils;
import application.InputParameters;
import application.PrecursorSet;

public class Sasita extends PrecursorFinder {

	OptimisationInterface ointerface;
	List<Compound> allCompounds = new LinkedList<Compound>();
	List<Reaction> allReactions = new LinkedList<Reaction>();
	//private int backupPartialSolutionsDelay = 300;
	/**
	 * The tolerance for consider that we have a positive production of the target.
	 */
	private double epsilon1 = 0.1;
	/**
	 * The tolerance for consider that we have a positive flux in a reaction i the MD-Model.
	 */
	private double epsilon2 = 0.1;
	/**
	 * The bigM value for the MILP models.
	 */
	private double bigM = 1000;
	private boolean duplicatingMachineryActive = false;
	private boolean steadyState = false;


	public Sasita(MetabolicNetwork network, boolean specialEmptySet) {
		super(network, specialEmptySet);
		allReactions.addAll(this.network.getReactions().values());
		Set<Compound> compounds = new HashSet<Compound>();
		for (Reaction r : allReactions) {
			compounds.addAll(r.getProduces().values());
			compounds.addAll(r.getSubstrates().values());
		}
		allCompounds.addAll(compounds);
		
	}

	public double getEpsilon1() {
		return epsilon1;
	}

	public void setEpsilon1(double epsilon1) {
		this.epsilon1 = epsilon1;
	}

	public double getEpsilon2() {
		return epsilon2;
	}

	public void setEpsilon2(double epsilon2) {
		this.epsilon2 = epsilon2;
	}

	public double getBigM() {
		return bigM;
	}

	public void setBigM(double bigM) {
		this.bigM = bigM;
	}

	public boolean isDuplicatingMachineryActive() {
		return duplicatingMachineryActive;
	}

	public void setDuplicatingMachineryActive(boolean duplicatingMachinery) {
		this.duplicatingMachineryActive = duplicatingMachinery;
	}

	public boolean isSteadyState() {
		return steadyState;
	}

	public void setSteadyState(boolean steadyState) {
		this.steadyState = steadyState;
	}

	// Method to find the set of minimal precursor sets that produce a target
	// compound
	// directly in the metabolic network
	public List<PrecursorSet> findPrecursorsInNetwork(List<Compound> targets) {
		allReactions.clear();
		allReactions.addAll(network.getReactions().values());
		// End of filtering empty reactions for the CDE network

		// Saving the processed network
		new MetabolicNetworkSBMLWriter("/tmp/preporcessedNetwork.xml")
				.write(this.getNetwork());

		System.out.println("After processing the network has "
				+ network.getReactions().size() + " reactions and "
				+ network.getCompounds().size() + " compounds.");

		List<Compound> prec = new LinkedList<Compound>();
		for (Compound c : this.network.getCompounds().values()) {
			if (c.isPrecursor()) {
				prec.add(c);
			}
		}
		System.out.println("Quantity of sources: " + prec.size());
		System.out
				.println("ID of all compounds considered as sources: " + prec);

		List<PrecursorSet> solutions = null;

		if (InputParameters.oneByOne) {
			// for each target, call
			// findPrecursorsInNetworkForTargetCompound(target)
			for (int j = 0; j < targets.size(); j++) {
				// starts computing the time to compute the solutions
				timeStart = System.currentTimeMillis();

				// Get the j-th target from targets and make it a target
				// and maps it to the restored network
				Compound target = getNetwork().getCompounds().get(
						targets.get(j).getId());
				target.setTarget(true);

				// writes down the reduced network for analysis
				try {
					new MetabolicNetworkSBMLWriter("preprocessedFor"
							+ target.getId() + ".xml").write(getNetwork());
				} catch (NullPointerException e) {
					e.printStackTrace();
					System.exit(-1);
				}
				
				System.out.println("Searching for precursors for targets: " + target);

				// Compute the precursors for this target
				solutions = findPrecursorsInNetworkForTarget(target,
						InputParameters.minimalityCheck);

				logln("Processing finished in "
						+ (System.currentTimeMillis() - timeStart)
						+ " ms.\n--------------------\n");

				if (solutions.size() == 0) {
					logln("No solution for: " + targets.get(j));
				}
				printSolutionsPerTarget(solutions, targets.get(j));
			}
		} else {
			// starts computing the time to compute the solutions
			Long start = System.currentTimeMillis();

			List<Reaction> networkReactions = new ArrayList<Reaction>();
			for (Reaction reac : this.network.getReactions().values()) {
				networkReactions.add(reac);
			}

			// Compute the precursors for all targets at once
			System.out.println("Searching for precursors for targets: " + targets);

			Compound target = createArtificialTargetCompound(targets);
			this.allReactions.clear();
			this.allReactions.addAll(this.network.getReactions().values());
			
			solutions = findPrecursorsInNetworkForTarget(target,
					InputParameters.minimalityCheck);

			if (solutions.size() == 0) {
				logln("No solution for: " + targets);
			}

			logln("Processing finished in "
					+ (System.currentTimeMillis() - start)
					+ " ms.\n--------------------\n");

			printSolutions(solutions, this.targets);
		}

		return solutions;

	}

	// Method to find the set of minimal precursor sets that produce a target
	// compound
	// directly in the metabolic network
	private List<PrecursorSet> findPrecursorsInNetworkForTarget(
			Compound target, boolean minimalityCheck) {
		List<Compound> sources = getMarkedSources(target);
		List<PrecursorSet> solutions = new LinkedList<PrecursorSet>();
		//solutions = loadPartialSolutions("DEBUG_temporary_solutions2.txt");
		//String partialSolutionsFilename = "partialSolutionsFor"
		//		+ target.getId() + ".txt";
        // This is not working in this version
		//solutions = loadPartialSolutions(partialSolutionsFilename);
		PrecursorSet lastSolution = null;
		List <PrecursorSet> lastSolutions = null;
		//int lastBackup = 0;
		// DEBUG
		//PrintWriter writer = null;
		//try {
			//writer = new PrintWriter("DEBUG_temporary_solutions.txt", "UTF-8");
		//} catch (Exception e){
		//	e.printStackTrace();
		//	System.exit(-9);
		//}
		this.ointerface = new CPLEXInterface(this.network);

		this.ointerface.startup();
		int numberOfSolutions = 0;
		int sizeOfBiggestsolution = 0;
		sizeOfBiggestsolution = 0;
		boolean searchForMultipleSolutions = true;
		boolean areThereMoreSolutions = true;
		while (areThereMoreSolutions) {
			PrecursorSet solution = null;
			List <PrecursorSet> multipleSolutions = null;

			if (!InputParameters.verbose){
				System.out.print("\rPrecursor sets found so far: " + numberOfSolutions + " ");
			}
			if (duplicatingMachineryActive) {
				searchForMultipleSolutions = false;
				solution = findNextDuplicatingMachineryMinimalPrecursor(
						sources, target, lastSolution);
				multipleSolutions = new LinkedList<PrecursorSet>();
				multipleSolutions.add(solution);
			} else if (steadyState) {
				searchForMultipleSolutions = false;
				solution = findNextSteadyStateMinimalPrecursor(sources, target,
						lastSolution);
				multipleSolutions = new LinkedList<PrecursorSet>();
				multipleSolutions.add(solution);
			} else {
				if (searchForMultipleSolutions){
					multipleSolutions = findNextsMinimalPrecursor(sources, target,lastSolutions);
					if (multipleSolutions == null){
						searchForMultipleSolutions = false;
						continue;
					}
				}
				else{
					solution = findNextMinimalPrecursor(sources, target,lastSolution);
					multipleSolutions = new LinkedList<PrecursorSet>();
					if (solution != null){
						multipleSolutions.add(solution);
					}
				}
			}
			Iterator<PrecursorSet> iteratorMultipleSolutions;
			if (multipleSolutions != null && multipleSolutions.size() > 0){
				lastSolutions = multipleSolutions;
				iteratorMultipleSolutions = multipleSolutions.iterator();
				while (iteratorMultipleSolutions.hasNext()) {
					solution = iteratorMultipleSolutions.next();
					if (solution != null && solution.getPrecursors().size() >= 0) {
						//writer.println(solution.getPrecursors().toString());
						//writer.flush();
						if (solution.getPrecursors().size() < sizeOfBiggestsolution){
							System.err.println("This solution can not be a solution, problem!");
							System.err.println(solution);
							System.err.println("Aborting!");
							System.exit(-1);
						}
						numberOfSolutions++;
						sizeOfBiggestsolution = Math.max(sizeOfBiggestsolution,solution.getPrecursors().size());
						solutions.add(solution);
						lastSolution = solution;
						logln("Solution found for " + target + ": "
								+ solution.getPrecursors());
						// This is not working in this version
						/*
						if (solutions.size() / backupPartialSolutionsDelay > lastBackup) {
							lastBackup = solutions.size() / backupPartialSolutionsDelay;
						savePartialSolutions(solutions, partialSolutionsFilename);
						}
						 */
					}
					else if (!searchForMultipleSolutions){ 
						areThereMoreSolutions = false;
						break;
					}
				} 
			}else if (!searchForMultipleSolutions) {
				areThereMoreSolutions = false;
			}
		}
		// This is not working in this version
		//cleanPartialSolutions(partialSolutionsFilename);
		this.ointerface.finish();
		//writer.close();
		return solutions;
	}

	private List<PrecursorSet> loadPartialSolutions(String filename) {
		List<PrecursorSet> solutions = new LinkedList<PrecursorSet>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(filename));
			String line = in.readLine();
			while (line != null) {
				line = line.replace("[", "");
				line = line.replace("]", "");
				StringTokenizer tokens = new StringTokenizer(line, ",");
				PrecursorSet p = new PrecursorSet();
				while (tokens.hasMoreTokens()) {
					String compoundId = tokens.nextToken();
					compoundId = compoundId.replace("<", "");
					compoundId = compoundId.replace(">", "");
					compoundId = compoundId.trim();
					Compound c = this.network.getCompounds().get(compoundId);
					p.addPrecursor(c);
				}
				solutions.add(p);
				line = in.readLine();
			}
			in.close();
		} catch (FileNotFoundException e) {
			logln("File with previous solutions not found, one will be created.");
		} catch (IOException e) {
			e.printStackTrace();
		}
		return solutions;
	}

	/*private void savePartialSolutions(List<PrecursorSet> solutions,
			String filename) {
		PrintWriter writer;
		try {
			writer = new PrintWriter(filename, "UTF-8");
			for (PrecursorSet p : solutions) {
				writer.println(p.getPrecursors().toString());
			}
			writer.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		}
	}

	private void cleanPartialSolutions(String filename) {
		try {
			File f = new File(filename);
			f.delete();
		} catch (Exception x) {
			System.err.println(x);
		}
	}*/

	private PrecursorSet findNextMinimalPrecursor(List<Compound> sources,
			Compound target, PrecursorSet lastSolution) {
		return this.ointerface.findNextMinimalPrecursor(sources, target, lastSolution, this.bigM, this.epsilon1);
	}
	
	private List<PrecursorSet> findNextsMinimalPrecursor(List<Compound> sources,
			Compound target, List<PrecursorSet> lastSolutions) {
		return this.ointerface.findNextsMinimalPrecursor(sources, target, lastSolutions, this.bigM, this.epsilon1);
	}
	private PrecursorSet findNextDuplicatingMachineryMinimalPrecursor(
			List<Compound> sources, Compound target, PrecursorSet lastSolution) {
		return this.ointerface.findNextDuplicatingMachineryMinimalPrecursor(sources, target, lastSolution, this.bigM, this.epsilon1, this.epsilon2);
	}

	private PrecursorSet findNextSteadyStateMinimalPrecursor(
			List<Compound> sources, Compound target, PrecursorSet lastSolution) {
		return this.ointerface.findNextSteadyStateMinimalPrecursor(sources, target, lastSolution, this.bigM, this.epsilon1);
	}

	public void printSolutions2Xml(List<PrecursorSet> setOfSets, String target) {

		PrintWriter writer;
		try {
			System.out.println("\nPrint solutions to " + target
					+ "_PS.xml");
			writer = new PrintWriter(target + "_PS.xml", "UTF-8");
			writer.println("<precursorSets>");
			int countingSets = 1;
			for (PrecursorSet set : setOfSets) {
				writer.println("\t<precursorSet id=\"" + countingSets++ + "\">");
				for (Compound source : set.getPrecursors()) {
					writer.println("\t\t<source id=\"" + source.getId()
							+ "\" name=\"" + source.getName()
							+ "\" compartment=\"" + source.getCompartment()
							+ "\" />");
				}
				writer.println("\t</precursorSet>");
			}
			writer.println("</precursorSets>");
			writer.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		}
	}

	private List<Compound> getMarkedSources(Compound target) {
		List<Compound> prec = new LinkedList<Compound>();
		for (Compound c : this.network.getCompounds().values()) {
			if (c.isPrecursor()) {
				prec.add(c);
			}
		}
		return prec;
	}

	protected Compound createArtificialTargetCompound(List<Compound> targets) {
		// Create an special reaction that takes as substrates the targets
		// passed
		// and that produces an special, TARGET, compound.
		Compound target = getNetwork()
				.addCompound("TARGET", "TARGET", "TARGET");
		Reaction special = getNetwork().addNewReaction(
				"SpecialReactionThatProducesTarget",
				"SpecialReactionThatProducesTarget", false);
		special.addProduct(target, 1.0);
		// special.addSubstrates(targets);
		for (Compound c : targets) {
			special.addSubstrate(c, 1.0);
		}

		// Create an other special reaction that takes as substrates the new
		// target created
		Compound targetMinimal = getNetwork().addCompound("TARGET_MINIMAL",
				"TARGET_MINIMAL", "TARGET_MINIMAL");
		// and that produces an special, TARGET MINIMAL, compound.
		Reaction specialMinimal = getNetwork().addNewReaction(
				"SpecialReactionThatProducesTargetMinimal",
				"SpecialReactionThatProducesTargetMinimal", false);
		specialMinimal.addProduct(targetMinimal, 1.0);
		specialMinimal.addSubstrate(target, 1.0);

		return targetMinimal;
	}

	public void printSolutions(List<PrecursorSet> setOfSets, List<Compound> tcs) {
		System.out.println();
		if (setOfSets == null) {
			System.out.println("No solution found.");
			return;
		}

		// first, order each solution
		for (PrecursorSet set : setOfSets) {
			Collections.sort(set.getPrecursors(), null);
		}

		// then order the set of solutions
		Collections.sort(setOfSets, null);

		int countingSets = 1;
		for (PrecursorSet set : setOfSets) {
			String separator = "";
			System.out.println(countingSets++ + ")");
			System.out.print("Precursors : { ");
			for (Compound compound : set.getPrecursors()) {
				String id = compound.getId();
				id = StringUtils.sbmlDecode(compound.getId());

				System.out.print(separator + id);
				separator = ", ";
			}
			System.out.println(" }");

			System.out.print("Supplementary bootstrap compounds : { ");
			separator = "";
			for (Compound compound : set.getBootstraps()) {
				String id = compound.getId();
				id = StringUtils.sbmlDecode(compound.getId());

				System.out.print(separator + id);
				separator = ", ";

			}
			System.out.println(" }");

			if (set.getReactions() != null && set.getReactions().size() > 0){
				System.out.print("Reactions : { ");
				separator = "";
				for (Reaction r : set.getReactions()) {
					String id = r.getId();
					System.out.print(separator + id);
					separator = ", ";

				}
				System.out.println(" }");
			}

			if (set.getCumulatedCompounds() != null && set.getCumulatedCompounds().size() > 0){
				System.out.print("Cumulated compounds : { ");
				separator = "";
				for (Compound c : set.getCumulatedCompounds()) {
					String id = c.getId();
					System.out.print(separator + id);
					separator = ", ";

				}
				System.out.println(" }");
			}

			if (set.getPrecursors().isEmpty() && set.getReactions() != null) {
				System.out.println("\n-----\n\nReaction in the cycle:\n");
				for (Reaction r : set.getReactions()) {
					System.out.println(r.getId());
				}
				System.out.println("\n\n----\n");
			}

		}
	}

	protected void printSolutionsPerTarget(List<PrecursorSet> setOfSets,
			Compound target) {
		System.out.println();
		if (setOfSets == null) {
			System.out.println("No solution found.");
			return;
		}

		// first, order each solution
		for (PrecursorSet set : setOfSets) {
			Collections.sort(set.getPrecursors(), null);
		}

		// then order the set of solutions
		Collections.sort(setOfSets, null);

		// output to xml file
		printSolutions2Xml(setOfSets, target.getId());

		int countingSets = 1;
		for (PrecursorSet set : setOfSets) {
			String separator = "";
			System.out.println(countingSets++ + ")");
			System.out.print("Precursors : { ");
			for (Compound compound : set.getPrecursors()) {
				String id = compound.getId();
				id = StringUtils.sbmlDecode(compound.getId());

				System.out.print(separator + id);
				separator = ", ";
			}
			System.out.println(" }");

			if (set.getStopCompounds().size() > 0) {
				System.out.print("Need solutions of target(s) : { ");
				separator = "";
				for (Compound compound : set.getStopCompounds()) {
					String id = compound.getId();
					id = StringUtils.sbmlDecode(compound.getId());

					System.out.print(separator + id);
					separator = ", ";
				}
				System.out.println(" }");
			}

			System.out.print("Supplementary bootstrap compounds : { ");
			separator = "";
			for (Compound compound : set.getBootstraps()) {
				String id = compound.getId();
				id = StringUtils.sbmlDecode(compound.getId());

				System.out.print(separator + id);
				separator = ", ";

			}
			System.out.println(" }");

			System.out.print("Reactions : { ");
			separator = "";
			for (Reaction r : set.getReactions()) {
				String id = r.getId();
				System.out.print(separator + id);
				separator = ", ";

			}
			System.out.println(" }");

			if (set.getPrecursors().isEmpty() && set.getReactions() != null) {
				System.out.println("\n-----\n\nReaction in the cycle:\n");
				for (Reaction r : set.getReactions()) {
					logln(r.getId());
				}
				System.out.println("\n\n----\n");
			}

		}

		Set<Reaction> reac = new HashSet<Reaction>();
		for (PrecursorSet set : setOfSets) {
			reac.addAll(set.getReactions());
		}

		System.out.println("\n\nUnion of reactions: " + reac.toString()
				+ "\n\n");

	}

	private void logln(String message) {
		if (InputParameters.verbose) {
			System.out.println("[SASITA] " + message);
		} else {
			System.out.print(".");
		}
	}

	public void checkMinimality(List<PrecursorSet> solutions) {
		logln("Checking for the minimality of the solutions...");
		int solutionNumber = 1;
		for (PrecursorSet k : solutions) {
			List<Compound> keyList = k.getPrecursors();
			for (int j = solutionNumber; j < solutions.size(); ++j) {
				List<Compound> testList = solutions.get(j).getPrecursors();
				if (keyList.containsAll(testList)) {
					logln("Solution " + solutionNumber + " is not minimal");
				} else if (testList.containsAll(keyList)) {
					logln("Solution " + (j + 1) + " is not minimal");
				}
			}
			++solutionNumber;
		}
		logln("Done!");
	}

}
