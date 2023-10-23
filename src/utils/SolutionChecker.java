package utils;

import ilog.concert.IloException;
import ilog.concert.IloLPMatrix;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplexModeler;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import metabolicNetwork.Compound;
import metabolicNetwork.MetabolicNetwork;
import metabolicNetwork.Reaction;
import application.InputParameters;

public class SolutionChecker {

	private MetabolicNetwork network;
	private List<Compound> targets;
	private double maxFlux;
	private double minTargetProd;
	private Compound artificialTarget;
	private IloCplex cplex = null;
	private List<String> cplexColsIds = null;
	private List<Compound> allCompounds;


	/**
	 * Constructor 
	 * @param network
	 * @param targets
	 * @throws Exception 
	 */
	public SolutionChecker(MetabolicNetwork network, List<Compound> targets) throws Exception {
		this.network = network;
		for (Compound c : targets){
			if (!network.getCompounds().keySet().contains(c.getId())){
				throw new Exception("Invalid target " + c);
			}
		}
		this.targets = targets;
		this.artificialTarget = createAndAddArtificialTargetCompound(targets);
		this.maxFlux = InputParameters.bigM;
		this.minTargetProd = InputParameters.epsilon1;
		this.allCompounds = new ArrayList<Compound>(network.getCompounds().values());
	}

	public MetabolicNetwork getNetwork() {
		return network;
	}

	public List<Compound> getTargets() {
		return targets;
	}

	public void checkSolutions(String fileName){
		List<List<Compound>> setOfFalseSolutions = new ArrayList<List<Compound>>();
		List<List<Compound>> setOfSolutions = parseSolutionFile(fileName);
		//System.out.println("Amount of solutions to check: " + setOfSolutions.size());
		int checkedSolutions = 0;
		for (List<Compound> solution : setOfSolutions){
			double prod = doFBA(solution, this.allCompounds, this.artificialTarget, this.maxFlux, this.minTargetProd);
			if (prod <= 0){
				System.err.println("False solution: " + solution);
				System.exit(-3);
				setOfFalseSolutions.add(solution);
			}
			System.out.print("\rChecking solutions (" + (++checkedSolutions)+ "/"+setOfSolutions.size() +") [");
			int i;
			for(i=0; i<=10*(checkedSolutions/(float)setOfSolutions.size());++i){
				System.out.print(".");
			}
			for(int j = i; j <= 10; ++j){
				System.out.print(" ");
			}
			System.out.print("]");
		}
		if (setOfFalseSolutions.size() > 0){
			System.out.println("\nFalse solutions:");
			for (List<Compound> solution : setOfFalseSolutions){
				System.out.println(solution);
			}
			System.out.println("Total of " + setOfFalseSolutions.size() + " false solutions:");


		}
		else{
			System.out.println("\nAll solutions are valid.");
		}
	}

	public boolean checkSolution(List<Compound> solution) {
		double prod = doFBA(solution, this.allCompounds, this.artificialTarget, this.maxFlux, this.minTargetProd);
		if (prod <= 0){
			return false;
		}
		return true;
	}

	/* Create an special reaction that takes as substrates the targets
	 * and that produces an special, TARGET, compound.
	 * If TARGET already exists, return it.
	 *  
	 *  @return TARGET compound
	 */
	private Compound createAndAddArtificialTargetCompound(List<Compound> targets) {
		if (!this.getNetwork().getCompounds().containsKey("TARGET")){
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
			return target;
		}
		return this.getNetwork().getCompounds().get("TARGET");
	}

	private List<List<Compound>> parseSolutionFile(String fileName) {
		List<List<Compound>> precursorSetSolutions = new ArrayList<List<Compound>>();
		try {
			String line = "";
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			line = br.readLine().trim();			
			List<Compound> ps = new ArrayList<Compound>();
			while(line != null){
				line = line.trim();
				if(line.startsWith("<source id")){ // create source compound
					String id = line.replaceAll("(.*id=\")(.*?)(\"\\s.*)", "$2");
					if(network.getCompounds().containsKey(id)){
						Compound c = network.getCompounds().get(id);
						ps.add(c);
					}
					else{
						System.out.println("ERROR: Can not find in the network the source " + id);
						System.exit(0);
					}

				}
				if(line.startsWith("</precursorSet>")){ // new precursorSet
					precursorSetSolutions.add(ps);
					ps = new ArrayList<Compound>();
				}
				line = br.readLine();
			}
			br.close();
		} catch (FileNotFoundException e) {
			System.err.println(e.getMessage());
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
		return precursorSetSolutions;
	}

	private double doFBA(List<Compound> sources, List<Compound> freeToAccumulate, Compound target, double inOrOutMaxAmount, double minTargetProd) {
		double solutionValue = 0;
		int nVars = this.network.getReactions().values().size() + this.network.getCompounds().values().size();

		if (this.cplex == null){
			IloCplex cplex = null;
			IloNumVar[] x = null;
			String[] reacIds = new String[nVars];
			try {
				cplex = new IloCplex();
				// cplex.setOut(null); // no output to console

				int nbReac = 0;
				for (Reaction r : this.network.getReactions().values()) {
					reacIds[nbReac++] = r.getId();
				}
				for (Compound c : this.network.getCompounds().values()) {
					reacIds[nbReac++] = c.getId() + "_producer";
				}
				
				this.cplexColsIds = Arrays.asList(reacIds);
				IloCplexModeler modeler = new IloCplexModeler();
				x = modeler.numVarArray(nVars, 0.0, Double.MAX_VALUE, reacIds);
				cplex.add(x);

				for (int i=0; i<nVars; ++i){
					cplex.addGe(x[i], 0.0, "Bound_Inf_"+i);
					cplex.addLe(x[i], inOrOutMaxAmount, "Bound_Sup_"+i);
				}
				IloLinearNumExpr fobj = modeler.linearNumExpr();

				for (Compound c : network.getCompounds().values()) {
					IloLinearNumExpr expr = modeler.linearNumExpr();
					if (sources.contains(c)) {
						// If it is a source, add a reaction that produces it from
						// nothing
						int delta = this.cplexColsIds.indexOf(c.getId()+"_producer");
						expr.addTerm(1.0, x[delta]);
					}
					else if (c.isBootstrap()){
						// If it is a bootstrap, add a reaction that produces it from
						// nothing
						int delta = this.cplexColsIds.indexOf(c.getId()+"_producer");
						expr.addTerm(1.0, x[delta]);
					}
					for (int i = 0; i < this.network.getReactions().values().size(); i++) {
						Reaction r = this.network.getReactions().get(reacIds[i]);
						if (r.getProduces().values().contains(c)) { // reaction j
							// produces
							// compound i
							expr.addTerm(r.getProductStochiometricValue(c), x[i]);
						} else if (r.getSubstrates().values().contains(c)) { // reaction
							// j
							// consumes
							// compound
							// i
							expr.addTerm(-r.getSubstrateStochiometricValue(c), x[i]);
						}
					}

					// Mv = 0 except if free to accumulate
					if ((freeToAccumulate != null && freeToAccumulate.contains(c))
							|| target.equals(c)) {
						if (target.equals(c)){
							cplex.addGe(expr, minTargetProd, c.getId()+"_tgt");
						}
						else{
							cplex.addGe(expr, 0.0, c.getId()+"_qss");
						}
					} else {
						cplex.addEq(expr, 0.0, c.getId()+"_stdstt");
					}
					//TODO: Seems redundant, checking, soon to be erased.
					//if (sources.contains(c)){
					//	cplex.addLe(expr, inOrOutMaxAmount, c.getId()+"_src");
					//}
					if (target.equals(c)) {
						fobj.add(expr);
					}
				}
				cplex.addMaximize(fobj,"objective");

			} catch (IloException e) {
				System.err.println("Concert exception caught: " + e);
			}
			try {
				cplex.exportModel("/tmp/m.lp");
				this.cplex = new IloCplex();
				this.cplex.importModel("/tmp/m.lp");
				this.cplex.setOut(null);
				this.cplex.setWarning(null);
				this.cplex.setParam(IloCplex.IntParam.MIPKappaStats, 2);
				this.cplex.setParam(IloCplex.IntParam.MIPKappaStats, 2);
				this.cplex.setParam(IloCplex.DoubleParam.EpInt, 1e-9);
				this.cplex.setParam(IloCplex.DoubleParam.EpRHS, 1e-9);
			} catch (IloException e) {
				e.printStackTrace();
			}
		}
		else{
			// CPLEX was already initialized, let's only change the constrains
			IloLPMatrix lp = (IloLPMatrix)cplex.LPMatrixIterator().next();
			for (Compound c : sources){
				int j = this.cplexColsIds.indexOf(c.getId()+"_producer");
				int i = 2*nVars + (j-this.network.getReactions().values().size());
				try {
					lp.setNZ(i, j, 1.0);
				} catch (IloException e) {
					e.printStackTrace();
				}
			}
		}

		try {
			// solve
			if (this.cplex.solve()) {
				solutionValue = this.cplex.getObjValue();
			}
			// Erase the producer reactions for the next iteration
			IloLPMatrix lp = (IloLPMatrix)cplex.LPMatrixIterator().next();
			for (Compound c : sources){
				int j = this.cplexColsIds.indexOf(c.getId()+"_producer");
				int i = 2*nVars + (j-this.network.getReactions().values().size());
				lp.setNZ(i, j, 0.0);
			}
		} catch (IloException e) {
			e.printStackTrace();
		}
		return solutionValue;
	}

}
