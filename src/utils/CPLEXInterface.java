package utils;

import ilog.concert.IloConstraint;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplex.Status;
import ilog.cplex.IloCplex.UnknownObjectException;
import ilog.cplex.IloCplexModeler;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;

import metabolicNetwork.Compound;
import metabolicNetwork.MetabolicNetwork;
import metabolicNetwork.Reaction;
import utils.SasitaCPLEXModelHolder.SasitaModelType;
import application.InputParameters;
import application.PrecursorSet;

public class CPLEXInterface implements OptimisationInterface {

	private SasitaCPLEXModelHolder modelHolder;
	SolutionChecker checker;
	private MetabolicNetwork network;
	private List<Reaction> allReactions = new LinkedList<Reaction>();
	private List<Compound> allCompounds = new LinkedList<Compound>();

	public CPLEXInterface (MetabolicNetwork network){
		this.network = network;
		allReactions.addAll(this.network.getReactions().values());
		Set<Compound> compounds = new HashSet<Compound>();
		for (Reaction r : allReactions) {
			compounds.addAll(r.getProduces().values());
			compounds.addAll(r.getSubstrates().values());
		}
		allCompounds.addAll(compounds);
		try {
			this.checker = new SolutionChecker(this.network, InputParameters.getTargetCompounds());
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(-7);
		}
		//try {
		//	IloCplex cplex = new IloCplex();
		//	cplex.importModel("/tmp/test.lp");
		//	cplex.solve();
		//	System.exit(0);
		//} catch (IloException e) {
			// TODO Auto-generated catch block
			//e.printStackTrace();
			//System.exit(0);
		//}
		
	}
	
	@Override
	public PrecursorSet findNextMinimalPrecursor(List<Compound> sources,
			Compound target, PrecursorSet lastSolution, double bigM, double epsilon1) {
		PrecursorSet solution = null;
		IloCplex cplex = null;
		IloNumVar[] x = null;
		IloIntVar[] Ind = null;
		
		if (this.modelHolder.isModelSet() == true
				&& this.modelHolder.getModelType() == SasitaModelType.NORMAL) {
			x = modelHolder.getVarX();
			Ind = modelHolder.getVarInd();
			cplex = modelHolder.getCplex();

			IloCplexModeler modeler = new IloCplexModeler();

			// Excluding last known solution
			try {
				if (lastSolution != null) {
					IloLinearNumExpr expr;
					expr = modeler.linearNumExpr();
					for (Compound c : lastSolution.getPrecursors()) {
						int ind = sources.indexOf(c);
						expr.addTerm(1.0, Ind[ind]);
					}
					cplex.addLe(expr,
							lastSolution.getPrecursors().size() - 1.0,
							"SolutionExclusion" + lastSolution.getPrecursors());
				}
			} catch (IloException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return null;
			}
		} else {
			try {
				cplex = new IloCplex();
				List<Compound> bootstrap = InputParameters.getBootstrapCompounds();

				String[] reacNames = new String[this.allReactions.size()
						+ sources.size() + bootstrap.size()];
				int nbReac = 0;
				for (Reaction r : this.allReactions) {
					reacNames[nbReac++] = r.getId();
				}
				for (Compound c : sources) {
					reacNames[nbReac++] = c.getId() + "_producer";
				}
				for (Compound c: bootstrap){
					reacNames[nbReac++] = c.getId() + "_bootstrap_producer";
				}
				IloCplexModeler modeler = new IloCplexModeler();

				x = modeler.numVarArray(
						this.allReactions.size() + sources.size() + bootstrap.size(), 0.0,
						Double.MAX_VALUE, reacNames);
				cplex.add(x);
				for (int i=0; i<this.allReactions.size() + sources.size() + bootstrap.size(); ++i){
					cplex.addRange(0.0,x[i], bigM, "Bounds_"+i);
				}
				Ind = modeler.intVarArray(sources.size(), 0, 1);

				for (int i = 0; i < sources.size(); ++i) {
					Ind[i].setName(sources.get(i).getId() + "_IND");
				}
				cplex.add(Ind);

				for (Compound c : network.getCompounds().values()) {
					IloLinearNumExpr expr = modeler.linearNumExpr();
					if (sources.contains(c)) {
						int delta = sources.indexOf(c);
						expr.addTerm(1.0, x[this.allReactions.size() + delta]);
						// IFF Constraint for the flux of artificial reactions and 
						// objective function variables
						/*IloLinearNumExpr intConstraitLHS = modeler
								.linearNumExpr();
						intConstraitLHS.addTerm(1.0, x[this.allReactions.size()
								+ delta]);
						IloLinearNumExpr intConstraitRHS = modeler
								.linearNumExpr();
						intConstraitRHS.addTerm(bigM, Ind[delta]);
						cplex.addLe(intConstraitLHS, intConstraitRHS, c.getId()
								+ "FluxIfUsed");

						intConstraitLHS = modeler.linearNumExpr();
						intConstraitLHS.addTerm(1.0, Ind[delta]);
						intConstraitRHS = modeler.linearNumExpr();
						intConstraitRHS.addTerm(1.0, x[this.allReactions.size()
								+ delta]);
						cplex.addLe(intConstraitLHS, intConstraitRHS, c.getId()
								+ "_NotUsedNotCounted");
						*/
						// Same constraint but a formulation with two IFThen constraints						
						//cplex.add(cplex.ifThen(cplex.eq(Ind[delta], 0.0),cplex.le(x[this.allReactions.size() + delta], 0.0)));
						//cplex.add(cplex.ifThen(cplex.le(x[this.allReactions.size() + delta], 0.0),cplex.eq(Ind[delta], 0.0)));
						
						// Same constraint but a formulation with the indicator constraint scheme	
						cplex.addEq(cplex.sum(Ind[delta], cplex.le(x[this.allReactions.size() + delta], 0.0)), 1.0, reacNames[this.allReactions.size() + delta]+"_IC");
					}
					else if (c.isBootstrap()){
						int bdelta = bootstrap.indexOf(c);
						expr.addTerm(1.0, x[this.allReactions.size() + sources.size() + bdelta]);
						IloLinearNumExpr intConstraitLHS = modeler.linearNumExpr();
						intConstraitLHS.addTerm(1.0, x[this.allReactions.size()
								+ sources.size() + bdelta]);
						cplex.addLe(intConstraitLHS, bigM, c.getId()
								+ "BootstrapFluxUpperBound");
					}
					for (int i = 0; i < this.allReactions.size(); i++) {
						Reaction r = this.allReactions.get(i);
						if (r.getProduces().values().contains(c)) { // reaction
																	// j
																	// produces
																	// compound
																	// i
							expr.addTerm(r.getProductStochiometricValue(c),
									x[i]);
						} else if (r.getSubstrates().values().contains(c)) { // reaction
																				// j
																				// consumes
																				// compound
																				// i
							expr.addTerm(-r.getSubstrateStochiometricValue(c),
									x[i]);
						}
					}

					if (target.equals(c)) {
						cplex.addGe(expr, epsilon1, c.getId());
					} else if (!c.isBootstrap()){
						cplex.addGe(expr, 0.0, c.getId());
					}
				}

				// Excluding last known solution
				if (lastSolution != null) {
					IloLinearNumExpr expr = modeler.linearNumExpr();
					for (Compound c : lastSolution.getPrecursors()) {
						int ind = sources.indexOf(c);
						expr.addTerm(1.0, Ind[ind]);
					}
					cplex.addLe(expr,
							lastSolution.getPrecursors().size() - 1.0,
							"SolutionExclusion" + lastSolution.getPrecursors());
				}
				/*else{
					List<PrecursorSet> solutions = loadPartialSolutions("/tmp/teste");
					int cc = 0;
					for (PrecursorSet sol : solutions){
						IloLinearNumExpr expr = modeler.linearNumExpr();
						for (Compound c : sol.getPrecursors()) {
							int ind = sources.indexOf(c);
							expr.addTerm(1.0, Ind[ind]);
						}
						cplex.addLe(expr,
								sol.getPrecursors().size() - 1.0,
								"SolutionExclusion" + sol);
						++cc;
					}
					
				}*/

				IloLinearNumExpr fobj = modeler.linearNumExpr();
				for (IloNumVar Ii : Ind) {
					fobj.addTerm(1.0, Ii);
				}
				cplex.addMinimize(fobj);

				//cplex = new IloCplex();
				//cplex.importModel("/tmp/BUG_MIPEnumerationModel.lp");
				
				if (!InputParameters.verbose) {
					cplex.setOut(null);
					cplex.setWarning(null);
				}
				cplex.setParam(IloCplex.DoubleParam.EpInt, 1e-9);
				cplex.setParam(IloCplex.DoubleParam.EpRHS, 1e-9);
				cplex.setParam(IloCplex.IntParam.Threads,
						InputParameters.nbThreads);
				cplex.setParam(IloCplex.BooleanParam.NumericalEmphasis, true);

			} catch (IloException e) {
				System.err.println("Concert exception caught: " + e);
				return null;
			}
		}
		
		modelHolder.setModel(cplex, x, Ind, SasitaModelType.NORMAL);

		try {
			cplex.exportModel("/tmp/MIPEnumerationModel.lp");
			cplex.setParam(IloCplex.DoubleParam.TiLim, 1e+75);

			// solve
			if (cplex.solve()) {
				logln("Solution status = " + cplex.getStatus());
				logln("Solution value = " + cplex.getObjValue());

				solution = new PrecursorSet();
				double[] is = cplex.getValues(modelHolder.getVarInd());
				for (int i = 0; i < is.length; ++i) {
					if (((int) (is[i] + 0.5)) == 1.0) {
						Compound s = sources.get(i);
						solution.addPrecursor(s);
					}
				}
				
				if (InputParameters.fbaCheck) {
					//double fbaResult = doFBA(solution.getPrecursors(),
					//		allCompounds, target, bigM, epsilon1);
					logln("Testing solution...");
					if (!checker.checkSolution(solution.getPrecursors())) {
						System.err
								.println("Error, the last solution found is not a true solution.");
						System.err
								.println("Choose a bigger epsilon or a smaller bigM.");
						return null;
					}
				}

				if (InputParameters.addReactionsToSolutions || InputParameters.addCumulatedCompoundsToSolutions) {
					Set<Reaction> reactionsInSolution = new HashSet<Reaction>();
					HashMap<String, Double> SV = new HashMap<String, Double>();
					double[] xs = cplex.getValues(modelHolder.getVarX());
					for (int j = 0; j < xs.length; ++j) {
						if (xs[j] > 0) {
							String rName = modelHolder.getVarX()[j].getName();
							for (Reaction r : allReactions) {
								if (r.getId().equals(rName)) {
									if (xs[j] > 1e-6) {
										reactionsInSolution.add(r);
										for (Compound p : r.getProduces().values()){
											double tmp = xs[j]*r.getProductStochiometricValue(p);
											tmp += SV.get(p.getId());
											SV.put(p.getId(), tmp);
										}
										for (Compound p : r.getSubstrates().values()){
											double tmp = -(xs[j]*r.getSubstrateStochiometricValue(p));
											tmp += SV.get(p.getId());
											SV.put(p.getId(), tmp);
										}
									}
								}
							}
						}
					}
					if (InputParameters.addReactionsToSolutions){
						solution.addReactions(reactionsInSolution);
					}
					if (InputParameters.addCumulatedCompoundsToSolutions){
						for(Compound c : allCompounds){
							if (SV.get(c.getId()) > 1e-6){
								solution.addCumullatedCompound(c);
							}
						}
					}
				}
				
				if (solution.getPrecursors().size() > this.modelHolder.getSizeOfBiggestSolution()){
					modelHolder.addSolutionSizeConstraint(solution.getPrecursors().size());
				}
			}
		} catch (IloException e) {
			System.err.print("Concert exception caught: ");
			e.printStackTrace();
			return null;
		}
		return solution;
	}

	@Override
	public PrecursorSet findNextDuplicatingMachineryMinimalPrecursor(
			List<Compound> sources, Compound target, PrecursorSet lastSolution, double bigM, double epsilon1, double epsilon2) {
		PrecursorSet solution = null;
		IloCplex cplex = null;
		IloNumVar[] x = null;
		IloIntVar[] Ind = null;

		if (this.modelHolder.isModelSet() == true
				&& this.modelHolder.getModelType() == SasitaModelType.DUPMACH) {
			x = modelHolder.getVarX();
			Ind = modelHolder.getVarInd();
			cplex = modelHolder.getCplex();

			IloCplexModeler modeler = new IloCplexModeler();

			// Excluding known solutions
			try {
				if (lastSolution != null) {
					IloLinearNumExpr expr;
					expr = modeler.linearNumExpr();
					for (Compound c : lastSolution.getPrecursors()) {
						int ind = sources.indexOf(c);
						expr.addTerm(1.0, Ind[ind]);
					}
					cplex.addLe(expr,
							lastSolution.getPrecursors().size() - 1.0,
							"SolutionExclusion" + lastSolution);
				}
			} catch (IloException e) {
				e.printStackTrace();
				return null;
			}
		} else {
			try {
				IloCplexModeler modeler = new IloCplexModeler();
				cplex = new IloCplex();
				List<Compound> bootstrap = InputParameters.getBootstrapCompounds();
				String[] reacNames = new String[this.allReactions.size()
						+ sources.size() + bootstrap.size()];
				int nbReac = 0;
				for (Reaction r : this.allReactions) {
					reacNames[nbReac++] = r.getId();
					for (Compound src : sources) {
						if (r.getProduces().containsKey(src.getId())) {
							r.getProduces().remove(src.getId());
						}
					}
				}
				for (Compound c : sources) {
					reacNames[nbReac++] = c.getId() + "_producer";
				}
				for (Compound c: bootstrap){
					reacNames[nbReac++] = c.getId() + "_bootstrap_producer";
				}
				x = modeler.numVarArray(
						this.allReactions.size() + sources.size() + bootstrap.size(), 0.0,
						Double.MAX_VALUE, reacNames);
				cplex.add(x);
				
				for (int i=0; i<this.allReactions.size() + sources.size() + bootstrap.size(); ++i){
					cplex.addRange(0.0,x[i], bigM, "Bounds_"+i);
				}

				Ind = modeler.intVarArray(sources.size(), 0, 1);

				for (int i = 0; i < sources.size(); ++i) {
					Ind[i].setName(sources.get(i).getId() + "_IND");
				}
				cplex.add(Ind);

				for (Compound c : network.getCompounds().values()) {
					int delta = sources.indexOf(c);
					IloLinearNumExpr expr = modeler.linearNumExpr();
					if (sources.contains(c)) {
						expr.addTerm(1.0, x[this.allReactions.size() + delta]);
						IloLinearNumExpr intConstraitLHS = modeler
								.linearNumExpr();
						intConstraitLHS.addTerm(1.0, x[this.allReactions.size()
								+ delta]);
						IloLinearNumExpr intConstraitRHS = modeler
								.linearNumExpr();
						intConstraitRHS.addTerm(bigM, Ind[delta]);
						cplex.addLe(intConstraitLHS, intConstraitRHS, c.getId()
								+ "FluxIfUsed");

						intConstraitLHS = modeler.linearNumExpr();
						intConstraitLHS.addTerm(1.0, Ind[delta]);
						intConstraitRHS = modeler.linearNumExpr();
						intConstraitRHS.addTerm(1.0, x[this.allReactions.size()
								+ delta]);
						cplex.addLe(intConstraitLHS, intConstraitRHS, c.getId()
								+ "_NotUsedNotCounted");

						// intConstraitLHS = modeler.linearNumExpr();
						// intConstraitLHS.addTerm(1.0, I[delta]);
						// intConstraitRHS = modeler.linearNumExpr();
						// intConstraitRHS.addTerm(1.0, Ind[delta]);
						// cplex.addGe(intConstraitLHS, intConstraitRHS,
						// c.getId()
						// + "Count");
					}
					else if (c.isBootstrap()){
						// TODO: Not working... correct this!
						int bdelta = bootstrap.indexOf(c);
						expr.addTerm(1.0, x[this.allReactions.size() + sources.size() + bdelta]);
						IloLinearNumExpr intConstraitLHS = modeler.linearNumExpr();
						intConstraitLHS.addTerm(1.0, x[this.allReactions.size()
								+ sources.size() + bdelta]);
						cplex.addLe(intConstraitLHS, bigM, c.getId()
								+ "BootstrapFluxUpperBound");
					}
					for (int i = 0; i < this.allReactions.size(); i++) {
						Reaction r = this.allReactions.get(i);
						if (r.getProduces().values().contains(c)) { // reaction
																	// j
																	// produces
																	// compound
																	// i
							expr.addTerm(r.getProductStochiometricValue(c),
									x[i]);
						} else if (r.getSubstrates().values().contains(c)) { // reaction
																				// j
																				// consumes
																				// compound
																				// i
							expr.addTerm(-r.getSubstrateStochiometricValue(c),
									x[i]);
						}
					}

					if (target.equals(c)) {
						cplex.addGe(expr, epsilon1, c.getId());
					} else {
						cplex.addGe(expr, 0.0, c.getId());
					}
				}

				// MD constraint
				// ----------------
				for (Compound c : network.getCompounds().values()) {
					IloLinearNumExpr C1 = modeler.linearNumExpr();
					IloLinearNumExpr C2 = modeler.linearNumExpr();
					if (!c.isPrecursor()) {
						if (c.isBootstrap()){
							int bdelta = bootstrap.indexOf(c);
							C1.addTerm(1.0, x[this.allReactions.size() + sources.size() + bdelta]);
						}
						for (int i = 0; i < this.allReactions.size(); i++) {
							Reaction r = this.allReactions.get(i);
							if (r.getProduces().values().contains(c)) { // reaction
																		// j
																		// produces
																		// compound
																		// i
								C1.addTerm(r.getProductStochiometricValue(c),
										x[i]);
							} else if (r.getSubstrates().values().contains(c)) { // reaction
																					// j
																					// consumes
																					// compound
																					// i
								C1.addTerm(
										-r.getSubstrateStochiometricValue(c),
										x[i]);
							}
						}
						IloConstraint CC1 = modeler.ge(C1, epsilon2, "KTMP_1_"
								+ c.getId());

						for (int i = 0; i < this.allReactions.size(); i++) {
							Reaction r = this.allReactions.get(i);
							if (r.getSubstrates().values().contains(c)) { // reaction
																			// j
																			// consumes
																			// compound
																			// i
								C2.addTerm(1.0, x[i]);
							}
						}
						IloConstraint CC2 = modeler.eq(C2, 0.0,
								"KTMP_2_" + c.getId());

						cplex.add(modeler.or(CC1, CC2,
								"MD_Constraint_" + c.toString()));
					}
				}

				// Excluding last known solution
				if (lastSolution != null) {
					IloLinearNumExpr expr = modeler.linearNumExpr();
					for (Compound c : lastSolution.getPrecursors()) {
						int ind = sources.indexOf(c);
						expr.addTerm(1.0, Ind[ind]);
					}
					cplex.addLe(expr,
							lastSolution.getPrecursors().size() - 1.0,
							"SolutionExclusion" + lastSolution);
				}

				IloLinearNumExpr fobj = modeler.linearNumExpr();
				for (IloIntVar Ii : Ind) {
					fobj.addTerm(1.0, Ii);
				}
				//cplex.addLe(fobj, 9, "SpecialConstraint");

				cplex.addMinimize(fobj);

				//cplex.exportModel("/tmp/MIPEnumerationModel.lp");
				cplex.setParam(IloCplex.IntParam.MIPKappaStats, 2);
				cplex.setParam(IloCplex.DoubleParam.EpInt, 1e-9);
				cplex.setParam(IloCplex.DoubleParam.EpRHS, 1e-9);

				if (!InputParameters.verbose) {
					cplex.setOut(null);
					cplex.setWarning(null);
				}
				cplex.setParam(IloCplex.IntParam.Threads,
						InputParameters.nbThreads);

				/*
				 * FileOutputStream cplexLogFile; try { cplexLogFile = new
				 * FileOutputStream("/tmp/cplex.log");
				 * cplex.setOut(cplexLogFile); // no output to console
				 * cplex.setWarning(cplexLogFile); } catch
				 * (FileNotFoundException e) { e.printStackTrace(); }
				 */
				// solve
			} catch (IloException e) {
				System.err.println(".Concert exception caught: " + e);
				return null;
			}
		}

		modelHolder.setModel(cplex, x, Ind, SasitaModelType.DUPMACH);

		try {
			if (cplex.solve()) {
				logln("Solution status = " + cplex.getStatus());
				logln("Solution value = " + cplex.getObjValue());

				solution = new PrecursorSet();
				double[] iss = cplex.getValues(Ind);
				for (int i = 0; i < iss.length; ++i) {
					iss[i] = (int) (iss[i] + 0.5);
					if (iss[i] == 1.0) {
						Compound s = sources.get(i);
						solution.addPrecursor(s);
					}
				}
			}
			if (InputParameters.addReactionsToSolutions) {
				Set<Reaction> reactionsInSolution = new HashSet<Reaction>();
				double[] xs = cplex.getValues(modelHolder.getVarX());
				for (int j = 0; j < xs.length; ++j) {
					if (xs[j] > 0) {
						String rName = modelHolder.getVarX()[j].getName();
						for (Reaction r : allReactions) {
							if (r.getId().equals(rName)) {
								if (xs[j] > 1e-6) {
									reactionsInSolution.add(r);
								}
							}
						}
					}
				}
				solution.addReactions(reactionsInSolution);
			}

			if (solution.getPrecursors().size() > this.modelHolder.getSizeOfBiggestSolution()){
				modelHolder.addSolutionSizeConstraint(solution.getPrecursors().size());
			}
			
		} catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
			return null;
		}

		return solution;
	}

	@Override
	public PrecursorSet findNextSteadyStateMinimalPrecursor(
			List<Compound> sources, Compound target, PrecursorSet lastSolution, double bigM, double epsilon1) {
		PrecursorSet solution = null;
		IloCplex cplex = null;
		IloNumVar[] x = null;
		IloIntVar[] Ind = null;
		
		if (this.modelHolder.isModelSet() == true
				&& this.modelHolder.getModelType() == SasitaModelType.STEADYSTATE) {
			x = modelHolder.getVarX();
			Ind = modelHolder.getVarInd();
			cplex = modelHolder.getCplex();

			IloCplexModeler modeler = new IloCplexModeler();
			// Excluding last known solution
			try {
				if (lastSolution != null) {
					IloLinearNumExpr expr;
					expr = modeler.linearNumExpr();
					for (Compound c : lastSolution.getPrecursors()) {
						int ind = sources.indexOf(c);
						expr.addTerm(1.0, Ind[ind]);
					}
					cplex.addLe(expr,
							lastSolution.getPrecursors().size() - 1.0,
							"SolutionExclusion" + lastSolution.getPrecursors());
				}
			} catch (IloException e) {
				e.printStackTrace();
				return null;
			}
		} else {
			try {
				cplex = new IloCplex();
				List<Compound> bootstrap = InputParameters.getBootstrapCompounds();
				List<Compound> boundary = new LinkedList<Compound>();
				for (Compound c : allCompounds) {
					if (c.isBoundary() || (!c.isBoundary() && sources.contains(c))){
						boundary.add(c);
					}
				}
				
				String[] reacNames = new String[this.allReactions.size()
				                				+ sources.size() + bootstrap.size() + boundary.size()];
				int nbReac = 0;
				for (Reaction r : this.allReactions) {
					reacNames[nbReac++] = r.getId();
				}
				for (Compound c : sources) {
					reacNames[nbReac++] = c.getId() + "_producer";
				}
				for (Compound c : bootstrap) {
					reacNames[nbReac++] = c.getId() + "_bootstrap";
				}
				for (Compound c : boundary) {
					reacNames[nbReac++] = c.getId() + "_boundary_EX";
				}
				IloCplexModeler modeler = new IloCplexModeler();
				x = modeler.numVarArray(
						this.allReactions.size() + sources.size() + bootstrap.size() + boundary.size(), 0.0,
						Double.MAX_VALUE, reacNames);
				cplex.add(x);
				
				for (int i=0; i<this.allReactions.size() + sources.size() + bootstrap.size() + boundary.size(); ++i){
						cplex.addRange(0.0,x[i], bigM, "Bounds_"+i);
				}
				
				Ind = modeler.intVarArray(sources.size(), 0, 1);

				for (int i = 0; i < sources.size(); ++i) {
					Ind[i].setName(sources.get(i).getId() + "_IND");
				}
				cplex.add(Ind);

				for (Compound c : network.getCompounds().values()) {
					IloLinearNumExpr expr = modeler.linearNumExpr();
					if (sources.contains(c)) {
						int delta = sources.indexOf(c);
						expr.addTerm(1.0, x[this.allReactions.size() + delta]);
						IloLinearNumExpr intConstraitLHS = modeler
								.linearNumExpr();
						intConstraitLHS.addTerm(1.0, x[this.allReactions.size()
								+ delta]);
						IloLinearNumExpr intConstraitRHS = modeler
								.linearNumExpr();
						intConstraitRHS.addTerm(bigM, Ind[delta]);
						cplex.addLe(intConstraitLHS, intConstraitRHS, c.getId()
								+ "FluxIfUsed");

						intConstraitLHS = modeler.linearNumExpr();
						intConstraitLHS.addTerm(1.0, Ind[delta]);
						intConstraitRHS = modeler.linearNumExpr();
						intConstraitRHS.addTerm(1.0, x[this.allReactions.size()
								+ delta]);
						cplex.addLe(intConstraitLHS, intConstraitRHS, c.getId()
								+ "_NotUsedNotCounted");
					}
					else if (c.isBootstrap()){
						int bdelta = bootstrap.indexOf(c);
						expr.addTerm(1.0, x[this.allReactions.size() + sources.size() + bdelta]);
						IloLinearNumExpr fluxConstraitLHS = modeler.linearNumExpr();
						fluxConstraitLHS.addTerm(1.0, x[this.allReactions.size()
								+ sources.size() + bdelta]);
						cplex.addLe(fluxConstraitLHS, bigM, c.getId()
								+ "BootstrapFluxUpperBound");
					}
					/*if (c.isBoundary()){
						int bdelta = boundary.indexOf(c);
						if (bdelta < 0){
							System.exit(-1);
						}
						expr.addTerm(-1.0, x[this.allReactions.size() + sources.size() + bootstrap.size() + bdelta]);
					}*/
					for (int i = 0; i < this.allReactions.size(); i++) {
						Reaction r = this.allReactions.get(i);
						if (r.getProduces().values().contains(c)) { // reaction
																	// j
																	// produces
																	// compound
																	// i
							expr.addTerm(r.getProductStochiometricValue(c),
									x[i]);
						} else if (r.getSubstrates().values().contains(c)) { // reaction
																				// j
																				// consumes
																				// compound
																				// i
							expr.addTerm(-r.getSubstrateStochiometricValue(c),
									x[i]);
						}
					}

					if (target.equals(c)) {
						cplex.addGe(expr, epsilon1, c.getId()+"_TARGET");
					} else if (sources.contains(c) || c.isBoundary()){
						cplex.addGe(expr, 0.0, c.getId()+"_ACC");
					} else {
						cplex.addEq(expr, 0.0, c.getId());	
					} 
				}
				// Excluding last known solution
				if (lastSolution != null) {
					IloLinearNumExpr expr = modeler.linearNumExpr();
					for (Compound c : lastSolution.getPrecursors()) {
						int ind = sources.indexOf(c);
						expr.addTerm(1.0, Ind[ind]);
					}
					cplex.addLe(expr,
							lastSolution.getPrecursors().size() - 1.0,
							"SolutionExclusion" + lastSolution.getPrecursors());
				}

				IloLinearNumExpr fobj = modeler.linearNumExpr();
				for (IloNumVar Ii : Ind) {
					fobj.addTerm(1.0, Ii);
				}
				cplex.addMinimize(fobj);

				if (!InputParameters.verbose) {
					cplex.setOut(null);
					cplex.setWarning(null);
				}
				cplex.setParam(IloCplex.DoubleParam.EpInt, 1e-9);
				cplex.setParam(IloCplex.DoubleParam.EpRHS, 1e-9);
				// cplex.exportModel("MIPEnumerationModel.lp");
				cplex.setParam(IloCplex.IntParam.Threads,
						InputParameters.nbThreads);
			} catch (IloException e) {
				System.err.println("Concert exception caught: " + e);
				return null;
			}
		}

		modelHolder.setModel(cplex, x, Ind, SasitaModelType.STEADYSTATE);

		try {
			// solve
			cplex.exportModel("/tmp/MIPEnumerationModel.lp");
			if (cplex.solve()) {
				logln("Solution status = " + cplex.getStatus());
				logln("Solution value = " + cplex.getObjValue());
				
				solution = new PrecursorSet();
				double[] is = cplex.getValues(modelHolder.getVarInd());
				for (int i = 0; i < is.length; ++i) {
					if (((int) (is[i] + 0.5)) == 1.0) {
						Compound s = sources.get(i);
						solution.addPrecursor(s);
					}
				}
				
				if (InputParameters.addReactionsToSolutions) {
					Set<Reaction> reactionsInSolution = new HashSet<Reaction>();
					double[] xs = cplex.getValues(modelHolder.getVarX());
					for (int j = 0; j < xs.length; ++j) {
						if (xs[j] > 0) {
							String rName = modelHolder.getVarX()[j].getName();
							for (Reaction r : allReactions) {
								if (r.getId().equals(rName)) {
									if (xs[j] > 1e-6) {
										reactionsInSolution.add(r);
									}
								}
							}
						}
					}
					solution.addReactions(reactionsInSolution);
				}
				
				if (solution.getPrecursors().size() > this.modelHolder.getSizeOfBiggestSolution()){
					modelHolder.addSolutionSizeConstraint(solution.getPrecursors().size());
				}
			}
		} catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
			return null;
		}
		return solution;
	}

	public List<PrecursorSet> findNextsMinimalPrecursor(List<Compound> sources,
			Compound target, List<PrecursorSet> lastSolutions, double bigM, double epsilon1) {
		List<PrecursorSet> solutions = null;
		IloCplex cplex = null;
		IloNumVar[] x = null;
		IloIntVar[] Ind = null;
		
		if (this.modelHolder.isModelSet() == true
				&& this.modelHolder.getModelType() == SasitaModelType.NORMAL) {
			x = modelHolder.getVarX();
			Ind = modelHolder.getVarInd();
			cplex = modelHolder.getCplex();

			IloCplexModeler modeler = new IloCplexModeler();

			// Excluding last known solutions
			try {
				if (lastSolutions != null) {
					for (PrecursorSet solution : lastSolutions){
						IloLinearNumExpr expr;
						expr = modeler.linearNumExpr();
						for (Compound c : solution.getPrecursors()) {
							int ind = sources.indexOf(c);
							expr.addTerm(1.0, Ind[ind]);
						}
						cplex.addLe(expr,
								solution.getPrecursors().size() - 1.0,
								"SolutionExclusion" + solution);
					}
				}
			} catch (IloException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return null;
			}
		} else {
			try {
				cplex = new IloCplex();
				List<Compound> bootstrap = InputParameters.getBootstrapCompounds();

				String[] reacNames = new String[this.allReactions.size()
						+ sources.size() + bootstrap.size()];
				int nbReac = 0;
				for (Reaction r : this.allReactions) {
					reacNames[nbReac++] = r.getId();
				}
				for (Compound c : sources) {
					reacNames[nbReac++] = c.getId() + "_producer";
				}
				for (Compound c: bootstrap){
					reacNames[nbReac++] = c.getId() + "_bootstrap_producer";
				}
				IloCplexModeler modeler = new IloCplexModeler();

				x = modeler.numVarArray(
						this.allReactions.size() + sources.size() + bootstrap.size(), 0.0,
						Double.MAX_VALUE, reacNames);
				cplex.add(x);
				for (int i=0; i<this.allReactions.size() + sources.size() + bootstrap.size(); ++i){
					cplex.addRange(0.0,x[i], bigM, "Bounds_"+i);
				}
				Ind = modeler.intVarArray(sources.size(), 0, 1);

				for (int i = 0; i < sources.size(); ++i) {
					Ind[i].setName(sources.get(i).getId() + "_IND");
				}
				cplex.add(Ind);

				for (Compound c : network.getCompounds().values()) {
					IloLinearNumExpr expr = modeler.linearNumExpr();
					if (sources.contains(c)) {
						int delta = sources.indexOf(c);
						expr.addTerm(1.0, x[this.allReactions.size() + delta]);
						// IFF Constraint for the flux of artificial reactions and 
						// objective function variables
						/*IloLinearNumExpr intConstraitLHS = modeler
								.linearNumExpr();
						intConstraitLHS.addTerm(1.0, x[this.allReactions.size()
								+ delta]);
						IloLinearNumExpr intConstraitRHS = modeler
								.linearNumExpr();
						intConstraitRHS.addTerm(bigM, Ind[delta]);
						cplex.addLe(intConstraitLHS, intConstraitRHS, c.getId()
								+ "FluxIfUsed");

						intConstraitLHS = modeler.linearNumExpr();
						intConstraitLHS.addTerm(1.0, Ind[delta]);
						intConstraitRHS = modeler.linearNumExpr();
						intConstraitRHS.addTerm(1.0, x[this.allReactions.size()
								+ delta]);
						cplex.addLe(intConstraitLHS, intConstraitRHS, c.getId()
								+ "_NotUsedNotCounted");
						*/
						// Same constraint but a formulation with two IFThen constraints						
						//cplex.add(cplex.ifThen(cplex.eq(Ind[delta], 0.0),cplex.le(x[this.allReactions.size() + delta], 0.0)));
						//cplex.add(cplex.ifThen(cplex.le(x[this.allReactions.size() + delta], 0.0),cplex.eq(Ind[delta], 0.0)));
						
						// Same constraint but a formulation with the indicator constraint scheme	
						cplex.addEq(cplex.sum(Ind[delta], cplex.le(x[this.allReactions.size() + delta], 0.0)), 1.0, reacNames[this.allReactions.size() + delta]+"_IC");
					}
					else if (c.isBootstrap()){
						int bdelta = bootstrap.indexOf(c);
						expr.addTerm(1.0, x[this.allReactions.size() + sources.size() + bdelta]);
						IloLinearNumExpr intConstraitLHS = modeler.linearNumExpr();
						intConstraitLHS.addTerm(1.0, x[this.allReactions.size()
								+ sources.size() + bdelta]);
						cplex.addLe(intConstraitLHS, bigM, c.getId()
								+ "BootstrapFluxUpperBound");
					}
					for (int i = 0; i < this.allReactions.size(); i++) {
						Reaction r = this.allReactions.get(i);
						if (r.getProduces().values().contains(c)) { // reaction
																	// j
																	// produces
																	// compound
																	// i
							expr.addTerm(r.getProductStochiometricValue(c),
									x[i]);
						} else if (r.getSubstrates().values().contains(c)) { // reaction
																				// j
																				// consumes
																				// compound
																				// i
							expr.addTerm(-r.getSubstrateStochiometricValue(c),
									x[i]);
						}
					}

					if (target.equals(c)) {
						cplex.addGe(expr, epsilon1, c.getId());
					} else if (!c.isBootstrap()){
						cplex.addGe(expr, 0.0, c.getId());
					}
				}

				// Excluding last known solutions
				if (lastSolutions != null) {
					for (PrecursorSet solution : lastSolutions){
						IloLinearNumExpr expr;
						expr = modeler.linearNumExpr();
						for (Compound c : solution.getPrecursors()) {
							int ind = sources.indexOf(c);
							expr.addTerm(1.0, Ind[ind]);
						}
						cplex.addLe(expr,
								solution.getPrecursors().size() - 1.0,
								"SolutionExclusion" + solution);
					}
				}
				/*else{
					List<PrecursorSet> solutions = loadPartialSolutions("/tmp/teste");
					int cc = 0;
					for (PrecursorSet sol : solutions){
						IloLinearNumExpr expr = modeler.linearNumExpr();
						for (Compound c : sol.getPrecursors()) {
							int ind = sources.indexOf(c);
							expr.addTerm(1.0, Ind[ind]);
						}
						cplex.addLe(expr,
								sol.getPrecursors().size() - 1.0,
								"SolutionExclusion" + sol);
						++cc;
					}
					
				}*/

				IloLinearNumExpr fobj = modeler.linearNumExpr();
				for (IloNumVar Ii : Ind) {
					fobj.addTerm(1.0, Ii);
				}
				cplex.addMinimize(fobj);

				//cplex = new IloCplex();
				//cplex.importModel("/tmp/BUG_MIPEnumerationModel.lp");
				
				if (!InputParameters.verbose) {
					cplex.setOut(null);
					cplex.setWarning(null);
				}
				cplex.setParam(IloCplex.DoubleParam.EpInt, 1e-9);
				cplex.setParam(IloCplex.DoubleParam.EpRHS, 1e-9);
				cplex.setParam(IloCplex.IntParam.Threads,
						InputParameters.nbThreads);
				cplex.setParam(IloCplex.BooleanParam.NumericalEmphasis, true);

			} catch (IloException e) {
				System.err.println("Concert exception caught: " + e);
				return null;
			}
		}
		
		modelHolder.setModel(cplex, x, Ind, SasitaModelType.NORMAL);

		try {
			cplex.exportModel("/tmp/MIPEnumerationModel.lp");
			// Set population parameters
		      cplex.setParam(IloCplex.DoubleParam.SolnPoolAGap, 0.5);
		      cplex.setParam(IloCplex.IntParam.SolnPoolIntensity, 2);
		      // Max amount of solution found by populate
		      cplex.setParam(IloCplex.IntParam.SolnPoolCapacity, 100);
		      cplex.setParam(IloCplex.IntParam.PopulateLim, 100);
		      cplex.setParam(IloCplex.DoubleParam.TiLim, InputParameters.TiLim);
		      cplex.setParam(IloCplex.IntParam.SolnPoolReplace, 1);			
			// solve
			if (cplex.populate()) {
				if (cplex.getStatus() != Status.Optimal){
					System.err.println("\nPopulate could not find a solution, changing strategy...");
					this.modelHolder.touchModelToResolve();
					return null;
				}
				logln("Solution status = " + cplex.getStatus());
				logln("Solution value = " + cplex.getObjValue());

				solutions = getSolutionsFromPopulate(cplex,sources,modelHolder.getSizeOfBiggestSolution());
				
				if (InputParameters.fbaCheck) {
					logln("Testing solutions...");
					List<PrecursorSet> falseSolutions = new ArrayList<PrecursorSet>();
					for (PrecursorSet solution : solutions){
						if (!checker.checkSolution(solution.getPrecursors())) {
							falseSolutions.add(solution);
						}
					}
					solutions.removeAll(falseSolutions);
					if (solutions.size() == 0){
						System.err
						.println("\nError, the last solutions found are not true solutions.");
						System.err
						.println("Trying again, now looking for one solution a time...");
						this.modelHolder.touchModelToResolve();
						return null;
					}
				}

				if (InputParameters.addReactionsToSolutions || InputParameters.addCumulatedCompoundsToSolutions) {
					for (PrecursorSet solution : solutions){
						Set<Reaction> reactionsInSolution = new HashSet<Reaction>();
						double[] xs = cplex.getValues(modelHolder.getVarX());
						HashMap<String, Double> SV = new HashMap<String, Double>();
						for (Compound c : allCompounds){
							SV.put(c.getId(), 0.0);
						}
						for (int j = 0; j < xs.length; ++j) {
							if (xs[j] > 0) {
								String rName = modelHolder.getVarX()[j].getName();
								for (Reaction r : allReactions) {
									if (r.getId().equals(rName)) {
										if (xs[j] > 1e-6) {
											reactionsInSolution.add(r);
											for (Compound p : r.getProduces().values()){
												double tmp = xs[j]*r.getProductStochiometricValue(p);
												tmp += SV.get(p.getId());
												SV.put(p.getId(), tmp);
											}
											for (Compound p : r.getSubstrates().values()){
												double tmp = -(xs[j]*r.getSubstrateStochiometricValue(p));
												tmp += SV.get(p.getId());
												SV.put(p.getId(), tmp);
											}
										}
									}
								}
							}
						}
						if(InputParameters.addReactionsToSolutions){
							solution.addReactions(reactionsInSolution);
						}
						if(InputParameters.addCumulatedCompoundsToSolutions){
							for(Compound c : allCompounds){
								if (SV.get(c.getId()) > 1e-6){
									solution.addCumullatedCompound(c);
								}
							}
						}
					}
				}
				for (PrecursorSet solution : solutions){
					if (solution.getPrecursors().size() > this.modelHolder.getSizeOfBiggestSolution()){
						modelHolder.addSolutionSizeConstraint(solution.getPrecursors().size());
					}
				}
			}
		} catch (IloException e) {
			System.err.print("Concert exception caught: ");
			e.printStackTrace();
			return null;
		}
		return solutions;
	}

	private List<PrecursorSet> getSolutionsFromPopulate(IloCplex cplex, List<Compound> sources, int minSizeAccepted) {
		List<PrecursorSet> solutions = new LinkedList<PrecursorSet>();
		int nsol = cplex.getSolnPoolNsolns();
		double smallObjFunction = Double.POSITIVE_INFINITY;
		List <Integer> realSolutions = new ArrayList<Integer>();
		try {
			for (int i = 0; i < nsol; i++) {
				if(cplex.getObjValue(i) > 0.2 && cplex.getObjValue(i) < 0.8 ){
					System.err.println("??? Aborting...");
					System.exit(-6);
				}
				int iObjF = (int) (cplex.getObjValue(i) + 0.5);
				if (iObjF >= minSizeAccepted && iObjF < smallObjFunction){
					smallObjFunction = iObjF;
					realSolutions.clear();
					realSolutions.add(i);
				} else if (iObjF >= minSizeAccepted && iObjF == smallObjFunction){
					realSolutions.add(i);
				}
			}
			for (int k : realSolutions){
				PrecursorSet solution = new PrecursorSet();
				double[] is;
				is = cplex.getValues(modelHolder.getVarInd(), k);
				for (int i = 0; i < is.length; ++i) {
					if (((int) (is[i] + 0.5)) == 1.0) {
						Compound s = sources.get(i);
						solution.addPrecursor(s);
					}
				}
				if (!solutions.contains(solution)){
					solutions.add(solution);
				}
			}
		} catch (UnknownObjectException e) {
			e.printStackTrace();
			System.exit(-2);
		} catch (IloException e) {
			e.printStackTrace();
			System.exit(-3);
		}
		return solutions;
	}

	private double doFBA(List<Compound> sources,
			List<Compound> freeToAccumulate, Compound target, double inOrOutMaxAmount, double minTargetProd) {
		double solutionValue = 0;
		IloCplex cplex = null;
		IloNumVar[] x = null;
		List<Compound> bootstrap = InputParameters.getBootstrapCompounds();
		String[] reacNames = new String[this.allReactions.size()
				+ sources.size() + bootstrap.size()];
		try {
			cplex = new IloCplex();
			// cplex.setOut(null); // no output to console

			int nbReac = 0;
			for (Reaction r : this.allReactions) {
				reacNames[nbReac++] = r.getId();
			}
			for (Compound c : sources) {
				reacNames[nbReac++] = c.getId() + "_producer";
			}
			IloCplexModeler modeler = new IloCplexModeler();
			x = modeler.numVarArray(this.allReactions.size() + sources.size() + bootstrap.size(),
					0.0, Double.MAX_VALUE, reacNames);
			cplex.add(x);

			for (int i=0; i<this.allReactions.size() + sources.size() + bootstrap.size(); ++i){
				cplex.addRange(0.0,x[i], inOrOutMaxAmount, "Bounds_"+i);
			}
			IloLinearNumExpr fobj = modeler.linearNumExpr();

			for (Compound c : network.getCompounds().values()) {
				IloLinearNumExpr expr = modeler.linearNumExpr();
				if (sources.contains(c)) {
					// If it is a source, add a reaction that produces it from
					// nothing
					int delta = sources.indexOf(c);
					expr.addTerm(1.0, x[this.allReactions.size() + delta]);
				}
				else if (c.isBootstrap()){
					// If it is a bootstrap, add a reaction that produces it from
					// nothing
					int delta = bootstrap.indexOf(c);
					expr.addTerm(1.0, x[this.allReactions.size() + sources.size() + delta]);
				}
				for (int i = 0; i < this.allReactions.size(); i++) {
					Reaction r = this.allReactions.get(i);
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
						cplex.addGe(expr, 0.0, c.getId());
					}
				} else {
					cplex.addEq(expr, 0.0, c.getId()+"_stdstt");
				}
				if (sources.contains(c)){
					cplex.addLe(expr, inOrOutMaxAmount, c.getId()+"_src");
				}
				if (target.equals(c)) {
					fobj.add(expr);
				}
			}
			cplex.addMaximize(fobj);
			if (!InputParameters.verbose) {
				cplex.setOut(null);
				cplex.setWarning(null);
			}
			cplex.exportModel("/tmp/fbaCheckModel.lp");

		} catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
			System.exit(-1);
		}
		
		try {
			// solve
			cplex.setParam(IloCplex.IntParam.MIPKappaStats, 2);
			cplex.setParam(IloCplex.DoubleParam.EpInt, 1e-9);
			cplex.setParam(IloCplex.DoubleParam.EpRHS, 1e-9);

			if (cplex.solve()) {
				logln("FBA test status = " + cplex.getStatus());
				solutionValue = cplex.getObjValue();
				logln("FBA solution value = " + solutionValue);
				if (solutionValue > 0) {
					logln("FBA OK");
				} else {
					logln("FBA Problem");
				}

				double[] xs = cplex.getValues(x);

				LinkedList<Reaction> reactionsInSolution = new LinkedList<Reaction>();
				for (int j = 0; j < xs.length; ++j) {
					if (xs[j] > 0) {
						String rName = x[j].getName();
						for (Reaction r : allReactions) {
							if (r.getId().equals(rName)) {
								if (xs[j] > 1e-6) {
									reactionsInSolution.add(r);
								}
							}
						}
					}
				}
			}
		} catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
			System.exit(-1);
		}
		return solutionValue;
	}

	private void logln(String message) {
		if (InputParameters.verbose) {
			System.out.println("[SASITA] " + message);
		} else {
			System.out.print(".");
		}
	}

	@Override
	public void startup() {
		this.modelHolder = new SasitaCPLEXModelHolder();
}

	@Override
	public void finish() {
		this.modelHolder.clearModel();		
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
}
