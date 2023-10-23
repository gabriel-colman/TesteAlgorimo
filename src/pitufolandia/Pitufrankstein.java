package pitufolandia;

import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplexModeler;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.TreeSet;

import metabolicNetwork.Compound;
import metabolicNetwork.MetabolicNetwork;
import metabolicNetwork.Reaction;
import pitufo.PrecursorFinder;
import utils.StringUtils;
import application.InputParameters;
import application.PrecursorSet;

/**
 * Implements the combinatorial version of the precursor set eumeration tool
 * described in the paper:
 * Eumeration of minimal stoichiometric precursor sets in metabolic networks
 * by Andrade, Wannagat et al.
 * 
 *
 */
public class Pitufrankstein extends PrecursorFinder {
	MetabolicNetwork manyToOneNetwork = new MetabolicNetwork();
	List<Reaction> allReactions = new LinkedList<Reaction>();
	LinkedList<Compound> allCompounds = new LinkedList<Compound>();
	/**
	 * The tolerance for consider that we have a positive production of the target.
	 */
	private double epsilon1 = 0.1;
	private double bigM = 1000.0;
	int maxK = 4;
	
	public Pitufrankstein(MetabolicNetwork network, boolean specialEmptySet) {
		super(network, specialEmptySet);
	}

	public double getEpsilon1() {
		return epsilon1;
	}

	public void setEpsilon1(double epsilon1) {
		this.epsilon1 = epsilon1;
	}

	/**
	 * @return the bigM
	 */
	public double getBigM() {
		return bigM;
	}

	/**
	 * @param bigM the bigM to set
	 */
	public void setBigM(double bigM) {
		this.bigM = bigM;
	}

	public List<PrecursorSet> findPrecursorsInNetworkForTarget(List<Compound> targets, boolean minimalityCheck)
	{
		Compound target = createArtificialTargetCompound(targets);
		List<Compound> sources = new ArrayList<Compound>();
		for (Compound c : this.network.getCompounds().values()){
			if (c.isPrecursor()){
				sources.add(c);
			}
		}
		System.out.println("Source compounds (" + sources.size() + "): " + sources);
		
		// Create artificial source compound for each original source compound.
		// This is important so that we don't get messed by internal reactions
		// that produce a source compound and so do't loose solutions
		List<Compound> artificialSources = new ArrayList<Compound>();
		for (Compound c : sources){
			c.setTopologicalPrecursor(false);
			c.setUserDefinedPrecursor(false);
			Compound nc = new Compound("___"+c.getId(), "___"+c.getId(), "");
			nc.setTopologicalPrecursor(true);
			nc.setUserDefinedPrecursor(true);
			Reaction nr = new Reaction("SR_"+c.getId(), "SR_"+c.getId(), false);
			nr.addSubstrate(nc,1.0);
			nr.addProduct(c,1.0);
			c.addProducedBy(nr);
			nc.addSubstrateOf(nr);
			this.network.getReactions().put(nr.getId(), nr);
			this.network.getCompounds().put(nc.getId(), nc);
			artificialSources.add(nc);
		}
		
		this.allReactions = new LinkedList<Reaction>(this.network.getReactions().values());
		this.allCompounds = new LinkedList<Compound>(this.network.getCompounds().values());
		
		// Generate the many-to-one network
		HashMap<String, Compound> ncompounds = new HashMap<String, Compound>();
		for (Compound c: allCompounds){
			Compound nc = new Compound(c.getId(),c.getName(),c.getCompartment());
			nc.setUserDefinedPrecursor(c.isPrecursor());
			ncompounds.put(nc.getId(), nc);
		}
		this.manyToOneNetwork.setCompounds(ncompounds);
		HashMap<String, Reaction> nreactions = new HashMap<String, Reaction>();
		for (Reaction r : allReactions){
			for (Compound t: r.getProduces().values()){
				Reaction nr = new Reaction(r.getId()+"_MTO_"+t.getId(),r.getName()+"_MTO_"+t.getName(),r.isReversible());
				for (String srcId : r.getSubstrates().keySet()){
					Compound e = ncompounds.get(srcId);
					e.getReactionsThatConsume(false).add(nr);
					nr.addSubstrate(e);
				}
				Compound e = ncompounds.get(t.getId());
				e.getProducedBy().add(nr);
				nr.addProduct(e);
				
				nreactions.put(nr.getId(), nr);
			}
		}
		for (Reaction r : nreactions.values()){
			String rId = r.getId();
			String rREVId;
			if (rId.substring(rId.length()-3, rId.length()).contains("_REV")){
				rREVId = rId.substring(0, rId.length()-3);
			}
			else{
				rREVId = rId+"_REV";
			}
			Reaction rREV = null;
			rREV = nreactions.get(rREVId);
			if (rREV != null){
				r.setReversible(true);
				rREV.setReversible(true);
				r.setReverse(rREV);
				rREV.setReverse(r);
			}
		}
		this.manyToOneNetwork.setReactions(nreactions);
		
		System.out.println("The many-to-one network has " + nreactions.size() + " reaction and " + ncompounds.size() + " compounds." );
		Compound targetManyToOne = this.manyToOneNetwork.getCompounds().get(target.getId());
		List <PrecursorSet> ps = visitCompound(targetManyToOne, null, new LinkedList<Compound>(), new TreeSet<Reaction>());
		// For each precursor set, replace the reactions from the many-to-one network 
		// by the corresponding original network reaction
		for (PrecursorSet p : ps){
			Set<Reaction> mtoreacs = p.getReactions();
			Set<Reaction> reacs = new HashSet<Reaction>();
			for(Reaction mtor : mtoreacs){
				String realID = mtor.getId();
				String[] tmp = realID.split("_MTO_");
				realID = tmp[0];
				Reaction realReac = this.network.getReactions().get(realID);
				reacs.add(realReac);
			}
			p.setReactions(reacs);
		}
		// Check feasibility
		List <PrecursorSet> solutions = new LinkedList<PrecursorSet>();
		List <PrecursorSet> combiner = new LinkedList<PrecursorSet>();
		for (PrecursorSet p : ps){
			double fluxMax = 0;
			fluxMax = doFBA(p.getPrecursors(), p.getBootstraps(), allCompounds, target);
			if (fluxMax > this.epsilon1){
				solutions.add(p);
			}
			else{
				combiner.add(p);
			}
		}
		System.out.println("Finished enumerating many-to-one tps");
		System.out.println("Combining " + combiner.size() + " tps solutions...");
		
		int k=2;
		while (combiner.size() > 1 && k <= maxK && !haveAllSolutions(solutions,artificialSources,target)){
			// Combine the unfeasible solutions up to size maxK and check feasibility
			List<PrecursorSet> powerSet = new LinkedList<PrecursorSet>(combiner);
			List <PrecursorSet> newcombiner = new LinkedList<PrecursorSet>();
			for (PrecursorSet o : combiner){
				for(PrecursorSet p : powerSet) {
					if (!p.getPrecursors().containsAll(o.getPrecursors())){
						PrecursorSet newp = new PrecursorSet(p);
						newp.add(o);
						double fluxMax = 0;
						fluxMax = doFBA(newp.getPrecursors(), newp.getBootstraps(), allCompounds, target);
						if (fluxMax > 0){
							if (!solutions.contains(newp)){
								solutions.add(newp);
							}
						}
						else{
							if (!newcombiner.contains(newp)){
								newcombiner.add(newp);
							}
						}
					}
				}
			}
			combiner = newcombiner;
			++k;
		}
		ps = null;
		
		// Remove non minimal precursor sets
		LinkedList<Integer> rl = new LinkedList<Integer>();
		for(int i = 0; i< solutions.size(); ++i){
			PrecursorSet p = solutions.get(i);
			p.setFlag(false);
			for(PrecursorSet pP : solutions){
				if (p != pP && !pP.isFlag() && precursorSetnIsSubsetOf(pP, p)){
					p.setFlag(true);
					rl.addFirst(i);
					break;
				}
			}
		}
		for(Integer i : rl){
			solutions.remove(i.intValue());
		}
		
		printSolutions(solutions);
		return solutions;
	}

	public List<PrecursorSet> visitCompound(Compound a, Reaction incommingReaction, List<Compound> visitedCompounds, Set<Reaction> visitedReactions){
//		if (incommingReaction != null) 
//			System.out.println("Visiting: " + a.getId() + "[Reaction: " + incommingReaction.getId() + "]");
		if (incommingReaction != null){
			visitedReactions.add(incommingReaction);
		}
		if (a.isPrecursor()){
			List<PrecursorSet> ps = new LinkedList<PrecursorSet>();
			PrecursorSet s = new PrecursorSet();
			s.addPrecursor(a);
			s.addReactions(visitedReactions);
			ps.add(s);
			return ps;
		}
		if (a.isBootstrap()){
			List<PrecursorSet> ps = new LinkedList<PrecursorSet>();
			PrecursorSet s = new PrecursorSet();
			s.addBootstrap(a);
			s.addReactions(visitedReactions);
			ps.add(s);
			return ps;
		}
		if (visitedCompounds.contains(a)){
			List<PrecursorSet> ps = new LinkedList<PrecursorSet>();
			PrecursorSet s = new PrecursorSet();
			s.addReactions(visitedReactions);
			s.addSubstrateAsStopCompound(a);
			ps.add(s);
			return ps;
		}

		visitedCompounds.add(a);

		Queue<Reaction> reactionsToAnalyze = new LinkedList<Reaction>();
		for(Reaction r : a.getProducedBy()){
			if (!visitedReactions.contains(r)){
				if((r.getReverseReaction() == null) || (r.getReverseReaction() != null && !visitedReactions.contains(r.getReverseReaction()))){
					reactionsToAnalyze.add(r);
				}
			}
		}
		

		List<PrecursorSet> allPSs = new LinkedList<PrecursorSet>();
		boolean allSolutionsFound = false;
		while( ! reactionsToAnalyze.isEmpty() && !allSolutionsFound){
			Reaction r = reactionsToAnalyze.poll();
			List<PrecursorSet> tps =  new LinkedList<PrecursorSet>();
			for(Compound c: r.getSubstrates().values()){
				List<Compound> newVisitedCompounds = new LinkedList<Compound>(visitedCompounds);
				List<Compound> compoundListToAdd = new LinkedList<Compound>(r.getSubstrates().values());
				compoundListToAdd.remove(c);
				newVisitedCompounds.addAll(compoundListToAdd);
				Set<Reaction> newVisitedReactions = new TreeSet<Reaction>(visitedReactions);
				List<PrecursorSet> lps = visitCompound(c,r, newVisitedCompounds, newVisitedReactions);
				tps = unionOfPSets(tps, lps);
			}
			for (PrecursorSet p : tps){
				if (!allPSs.contains(p))
					allPSs.add(p);
			}
		}

		// Minimalise the solution list
		List<PrecursorSet> reducedL = new ArrayList<PrecursorSet>(allPSs);
		for(int i = 0; i < allPSs.size(); i++)
		{
			for(int j = 0; j < allPSs.size(); j++)
			{
				if(i == j)
					continue;
				if( !allPSs.get(j).isFlag() && factoryIsSubsetOf(allPSs.get(j), allPSs.get(i)) )
				{
					reducedL.remove(allPSs.get(i));
					allPSs.get(i).setFlag( true );
					break;
				}
			}
		}
		allPSs = null;

		return reducedL;
	}

	protected Compound createArtificialTargetCompound(List<Compound> targets) {
		// Create an special reaction that takes as substrates the targets passed 
		// and that produces an special, TARGET, compound.
		Compound target = getNetwork().addCompound("TARGET", "TARGET", "TARGET");
		Reaction special = getNetwork().addNewReaction("SpecialReactionThatProducesTarget", "SpecialReactionThatProducesTarget", false);
		special.addProduct(target,1.0);
		for (Compound t: targets){
			special.addSubstrate(t,1.0);
		}
		
		// Create an other special reaction that takes as substrates the new target created 
		Compound targetMinimal = getNetwork().addCompound("TARGET_MINIMAL", "TARGET_MINIMAL", "TARGET_MINIMAL");
		// and that produces an special, TARGET MINIMAL, compound.
		Reaction specialMinimal = getNetwork().addNewReaction("SpecialReactionThatProducesTargetMinimal", "SpecialReactionThatProducesTargetMinimal", false);
		specialMinimal.addProduct(targetMinimal,1.0);
		specialMinimal.addSubstrate(target,1.0);		
		
		return targetMinimal;
	}
	
	static private boolean factoryIsSubsetOf(PrecursorSet s1, PrecursorSet s2)
	{
		HashSet<Compound> hash = new HashSet<Compound>();
		for(int i = 0; i < s2.getPrecursors().size(); i++)
			hash.add(s2.getPrecursors().get(i));
		
		boolean isAPrecusorSubset = hash.containsAll(s1.getPrecursors());
		
		if (isAPrecusorSubset){
			return s2.getReactions().contains(s1.getReactions());
		}
		
		return false;
	}	
	
	static private boolean precursorSetnIsSubsetOf(PrecursorSet s1, PrecursorSet s2)
	{
		HashSet<Compound> hash = new HashSet<Compound>();
		for(int i = 0; i < s2.getPrecursors().size(); i++)
			hash.add(s2.getPrecursors().get(i));
		
		boolean isAPrecusorSubset = hash.containsAll(s1.getPrecursors());
		
		return isAPrecusorSubset;
	}	
	
	/**
	 * Make the union of all PrecursorSets in lps with all Precursor Sets
	 * in tps.
	 */
	private List<PrecursorSet> unionOfPSets(List<PrecursorSet> tps, List<PrecursorSet> lps) {
		if (lps.size() > 0 ){
			if (tps.size() == 0){
				tps.addAll(lps);
				return tps;
			}
			else{
				List<PrecursorSet> newTps = new LinkedList<PrecursorSet>();
				for (PrecursorSet ps : lps){
					for (PrecursorSet x : tps){
						PrecursorSet newx = new PrecursorSet(x);
						newx.add(ps);
						boolean isRepeatedNewPrecursorSet = false;
						for(PrecursorSet s : newTps){
							if(s.haveSameSources(newx)){
								isRepeatedNewPrecursorSet = true;
								break;
							}
						}
						if(!isRepeatedNewPrecursorSet){
							newTps.add(newx);
						}
					}
				}
				return newTps;
			}
		}
		return tps;
	}
	private boolean haveAllSolutions(List<PrecursorSet> solutions, List<Compound> sources, Compound target) {
		boolean allSolutions = true;
		IloCplex cplex = null;
		IloNumVar[] x = null;
		IloIntVar[] Ind = null;
		String[] reacNames = new String[this.allReactions.size()
				+ sources.size()];
		
		try {
			cplex = new IloCplex();
	
			int nbReac = 0;
			for (Reaction r : this.allReactions) {
				reacNames[nbReac++] = r.getId();
			}
			for (Compound c : sources) {
				reacNames[nbReac++] = c.getId() + "_producer";
			}
			IloCplexModeler modeler = new IloCplexModeler();
			x = modeler.numVarArray(
					this.allReactions.size() + sources.size(), 0.0,
					Double.MAX_VALUE, reacNames);
			cplex.add(x);
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
					intConstraitRHS.addTerm(1000, Ind[delta]);
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
					cplex.addGe(expr, 0.1, c.getId());
				} else if (!c.isBootstrap()){
					cplex.addGe(expr, 0.0, c.getId());
				}
			}
	
			// Excluding last known solution
			for (PrecursorSet s : solutions){
				IloLinearNumExpr expr = modeler.linearNumExpr();
				for (Compound c : s.getPrecursors()) {
					int ind = sources.indexOf(c);
					expr.addTerm(1.0, Ind[ind]);
				}
				cplex.addLe(expr,
						s.getPrecursors().size() - 1.0,
						"SolutionExclusion" + s);
			}
			//IloLinearNumExpr fobj = modeler.linearNumExpr();
			//for (IloNumVar Ii : Ind) {
			//	fobj.addTerm(1.0, Ii);
			//}
			//cplex.addMinimize(fobj);
	
			if (!InputParameters.verbose) {
				cplex.setOut(null);
				cplex.setWarning(null);
			}
			cplex.setParam(IloCplex.DoubleParam.EpInt, 1e-9);
			cplex.setParam(IloCplex.DoubleParam.EpRHS, 1e-9);
			cplex.setParam(IloCplex.IntParam.Threads,
					InputParameters.nbThreads);
		} catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
		}
		
		
		try {
			if (cplex.solve()) {
				allSolutions = false;
			}
		} catch (IloException e) {
			System.out.println("Exception thronw, no solution found");
		}
		
		return allSolutions;
	}

	private void printSolutions(List<PrecursorSet> setOfSets) {
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

			System.out.print("Reactions (" + set.getReactions().size() + ") : { ");
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
					System.out.println(r.getId());
				}
				System.out.println("\n\n----\n");
			}

		}
	}

	private double doFBA(List<Compound> sources, List<Compound> bootstrap,
			List<Compound> freeToAccumulate, Compound target) {
		double solutionValue = 0;
		double inOrOutMaxAmount = this.bigM;
		IloCplex cplex = null;
		IloNumVar[] x = null;
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
					0.0, inOrOutMaxAmount, reacNames);
			cplex.add(x);
			
			for (int i=0; i<this.allReactions.size() + sources.size() + bootstrap.size(); ++i){
				cplex.addGe(x[i], 0.0, "Bound_Inf_"+i);
				cplex.addLe(x[i], inOrOutMaxAmount, "Bound_Sup_"+i);
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
					cplex.addGe(expr, 0.0, c.getId());
					cplex.addLe(expr, inOrOutMaxAmount, c.getId());
				} else {
					cplex.addEq(expr, 0.0, c.getId());
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
		} catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
		}

		try {
			// solve
			if (cplex.solve()) {
				solutionValue = cplex.getObjValue();
			}
		} catch (IloException e) {
			System.err.println("Concert exception caught: " + e);
		}
		return solutionValue;
	}
}

