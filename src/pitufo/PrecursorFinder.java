package pitufo;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import metabolicNetwork.Compound;
import metabolicNetwork.MetabolicNetwork;
import metabolicNetwork.Reaction;
import utils.MetabolicNetworkSBMLWriter;
import utils.StringUtils;
import application.InputParameters;
import application.PrecursorSet;

public class PrecursorFinder {

	boolean _STRICTED_GREATER_THAN = true;
	boolean _GREATER_THAN_OR_EQUAL = false;
	protected boolean specialEmptySet = false;		
	
	protected int numRecursions = 0;
	protected long timeStart = 0;

	protected MetabolicNetwork backup;
	protected MetabolicNetwork network;
	
	protected List<Compound> inputs = new ArrayList<Compound>();
	protected List<Compound> bootstraps = new ArrayList<Compound>();
	protected List<Compound> targets = new ArrayList<Compound>();
	protected List<Compound> userDefinedPrecursors = new ArrayList<Compound>();
	
	public PrecursorFinder(MetabolicNetwork network, boolean specialEmptySet)
	{
		this.specialEmptySet = specialEmptySet;
		defineWorkingNetwork(network);
	}	
	
	protected void restoreNetwork() {		
		setNetwork(backup.hardCopy());
		// update references
		updateReferences();
	}
	
	public void defineWorkingNetwork(MetabolicNetwork network) {
		this.network = network.hardCopy();
		this.backup = network.hardCopy();
		updateReferences();		
	}
	
	protected void updateReferences() {
		// update references
		targets.clear();
		for(Compound c : InputParameters.getTargetCompounds()) {
			targets.add( getNetwork().getCompounds().get(c.getId()));
		}
		
		inputs.clear();
		for(Compound c : InputParameters.getInputCompounds()) {
			inputs.add( getNetwork().getCompounds().get(c.getId()));
		}

		bootstraps.clear();
		for(Compound c : InputParameters.getBootstrapCompounds()) {
			bootstraps.add( getNetwork().getCompounds().get(c.getId()));
		}
		
		userDefinedPrecursors.clear();
		for(Compound c : InputParameters.getUserDefinedPrecursors()) {
			userDefinedPrecursors.add( getNetwork().getCompounds().get(c.getId()));
		}		
	}
	
	public MetabolicNetwork getNetwork() {
		return network;
	}

	public void setNetwork(MetabolicNetwork network) {
		this.network = network;
	}
	
	// Method to find the set of minimal precursor sets that produce a target compound
	// directly in the metabolic network
	public List<PrecursorSet> findPrecursorsInNetwork(List<Compound> targets)
	{
		System.out.println("\n\nSearching for precursors IN NETWORK\n");
		
		List<PrecursorSet> solutions = null;		
		
		if( InputParameters.oneByOne )
		{
			List<PrecursorSet> solution = null;
			List<Compound> individualTarget = new ArrayList<Compound>();
			// for each target, call findPrecursorsInNetworkForTarget(target)
			for(int j=0; j < targets.size(); j++)
			{
				// starts computing the time to compute the solutions
				timeStart = System.currentTimeMillis();

				// Get the j-th target from targets and make it a target
				// and maps it to the restored network
				Compound target = getNetwork().getCompounds().get(targets.get(j).getId());
				target.setTarget(true);

				// first decide which are the topological precursors
				network.extendTopologicalPrecursors( InputParameters.precursorIfProducedOnlyByReversible );
				
				// removes the "new" topological precursors obtained after reaction removals
				network.removeTopologicalSourcesThatAreNotPrecursors();
				System.out.println("\nThe reduced network after removing sources that are not precursors has "+network.getReactions().size()+" reactions and "+network.getCompounds().size()+" compounds.");								
				System.out.println("Num Precursors: "+getNetwork().getNumPrecursors());
				
				// writes down the reduced network for analysis
				new MetabolicNetworkSBMLWriter("preprocessedFor"+target.getId()+".xml").write(getNetwork());
				
				// Prepare a set of targets containing only the j-th target compound
				individualTarget.clear();
				individualTarget.add(target);
				
				// Compute the precursors for this target
				System.out.println("Searching for precursors for targets: "+individualTarget);
				solutions = findPrecursorsInNetworkForTarget( individualTarget, InputParameters.minimalityCheck);

				//printSolutions(solutions, this.targets);
				System.out.println("Processing finished in "+(System.currentTimeMillis()-timeStart)+" ms.\n--------------------\n");
				System.out.println("Number of recursions: "+numRecursions);

				// Restores a copy of the network for computing the results for the current target
				restoreNetwork();
			}
			return solution;
		}
		else
		{
			// starts computing the time to compute the solutions			
			Long start = System.currentTimeMillis();

			// Make all compounds in the targets set to be marked as targets
			for(int j=0; j < targets.size(); j++)
				targets.get(j).setTarget(true);
			
			// Compute the precursors for all targets at once
			System.out.println("Searching for precursors for targets: "+targets);
			List<PrecursorSet> solution = findPrecursorsInNetworkForTarget(targets, InputParameters.minimalityCheck);
			
			System.out.println("Processing finished in "+(System.currentTimeMillis()-start)+" ms.\n--------------------\n");
			return solution;
		}
	}	
	
	protected Compound createArtificialTargetCompound(List<Compound> targets) {
		// Create an special reaction that takes as substrates the targets passed 
		// and that produces an special, TARGET, compound.
		Compound target = getNetwork().addCompound("TARGET", "TARGET", "TARGET");
		Reaction special = getNetwork().addNewReaction("SpecialReactionThatProducesTarget", "SpecialReactionThatProducesTarget", false);
		special.addProduct(target);
		special.addSubstrates(targets);
		
		// Create an other special reaction that takes as substrates the new target created 
		Compound targetMinimal = getNetwork().addCompound("TARGET_MINIMAL", "TARGET_MINIMAL", "TARGET_MINIMAL");
		// and that produces an special, TARGET MINIMAL, compound.
		Reaction specialMinimal = getNetwork().addNewReaction("SpecialReactionThatProducesTargetMinimal", "SpecialReactionThatProducesTargetMinimal", false);
		specialMinimal.addProduct(targetMinimal);
		specialMinimal.addSubstrate(target);		
		
		return targetMinimal;
	}
	
	public List<PrecursorSet> findPrecursorsInNetworkForTarget(List<Compound> targets, boolean minimalityCheck)
	{
		return null;
	}
	
	public boolean ReactionIsSupersetOfSomeReaction(Reaction r, List<Reaction> reactionsAnalyzed, HashMap<String, Compound> Hk)
	{
		for(int i = 0; i < reactionsAnalyzed.size(); i++)
		{
			if( ReactionIsSupersetOf(r, reactionsAnalyzed.get(i), Hk, _STRICTED_GREATER_THAN) )
				return true;
		}
		return false;
	}
	
	public boolean ReactionIsSupersetOf(Reaction r1, Reaction r2, HashMap<String, Compound> Hk, boolean strictedGreater)
	{
		HashMap<String, Compound> substratesR1 = new HashMap<String, Compound>(r1.getSubstrates());
		HashMap<String, Compound> substratesR2 = new HashMap<String, Compound>(r2.getSubstrates());
		
		if( Hk != null)
		{
			for(Compound c: Hk.values()) {
				substratesR1.remove(c.getId());
				substratesR2.remove(c.getId());
			}
		}
			
		for(Compound inR2: substratesR2.values())
		{
			if( !substratesR1.containsKey( inR2.getId() ) )
				return false;
		}
		// check if they are the same...
		return !strictedGreater || (substratesR1.size() > substratesR2.size());
	}
	
/*	public List<Reaction> findMinimalReactionsThatProduce(Compound a, HashMap<String, Compound> Hk)
	{
		numIterations++;
		List<Reaction> result = new ArrayList<Reaction>();
		List<Reaction> reactionsThatProduce = a.getReactionsThatProduce(false);
		boolean isMinimal;
		for(int i = 0; i < reactionsThatProduce.size(); i++)
		{
			isMinimal = true;
			for(int j = 0; j < reactionsThatProduce.size(); j++)
			{
				if( i != j )
				{
					numComparisons++;
					if( ReactionIsSupersetOf(reactionsThatProduce.get(i), reactionsThatProduce.get(j), Hk, _GREATER_THAN_OR_EQUAL) )
					{
						if( reactionsThatProduce.get(i).getSubstrates().size() > reactionsThatProduce.get(j).getSubstrates().size() || i > j )
						{
							isMinimal = false;
							break;
						}
					}
				}
			}
			if( isMinimal )
				result.add(reactionsThatProduce.get(i));
		}
		return result;
	}*/
	
	public List<Reaction> findMinimalReactionsThatProduce(Compound a, HashMap<String, Compound> A)
	{
		List<Reaction> reactionsThatProduce = a.getReactionsThatProduce(false);
		if( reactionsThatProduce.size() <= 1 )
			return reactionsThatProduce;
		
		// result starts with the first reaction
		List<Reaction> result = new ArrayList<Reaction>();
		result.add(reactionsThatProduce.get(0));
		// and then try to add other reactions to the set
		for(int i = 1; i < reactionsThatProduce.size(); i++)
		{
			boolean isMinimal = false; boolean markedAsNonMinimal = false;
			for(int j = result.size()-1; j >= 0; j--)
			{
				if( ReactionIsSupersetOf(result.get(j), reactionsThatProduce.get(i), A, _STRICTED_GREATER_THAN) )
				{
					isMinimal = true;
					result.remove(j);
				}
				
				if( !isMinimal && ReactionIsSupersetOf(reactionsThatProduce.get(i), result.get(j), A, _GREATER_THAN_OR_EQUAL) ) {
					markedAsNonMinimal = true;
					break;
				}
			
			}
			if( isMinimal || !markedAsNonMinimal )
				result.add(reactionsThatProduce.get(i));
		}
		return result;
	}	
	
	public boolean setContainsReaction(List<Reaction> set, Reaction reaction, boolean checkReversible)
	{
		if( checkReversible )
		{
			if( reaction.isReversible() )
				return set.contains(reaction) || set.contains(reaction.getReverseReaction());
			else
				return set.contains(reaction);
		}
		else
			return set.contains(reaction);
	}
	
	protected Compound choosePivot(Reaction r, HashMap<String, Compound> A) {
		for(Compound c: r.getSubstrates().values()) {
			if( !c.isPrecursor() && !A.containsKey(c.getId()))
				return c;
		}
		return null;
	}
	
	public Reaction findNewMinimalReactionThatProduces(Compound target, HashMap<String, Compound> A)
	{
		List<Reaction> reactionsThatProduce = target.getReactionsThatProduce(false);
	
		for(int i = 0; i < reactionsThatProduce.size(); i++)
		{
			Compound pivot = choosePivot(reactionsThatProduce.get(i), A);
			if( pivot == null )
				continue;
			
			// The reaction "i" is the candidate to be the next new minimal reaction to be analyzed
			// Let's check if it is not contained in the rest of the reactions
			boolean isMinimal = true;
			int k = i;
			for(int j = k+1; j < reactionsThatProduce.size(); j++)
			{
				// If the candidate "i" is a superset of the reaction "j", then we change the candidate to "j" 
				if( ReactionIsSupersetOf(reactionsThatProduce.get(i), reactionsThatProduce.get(j), A, _STRICTED_GREATER_THAN) ) {
					Compound pivotJ = choosePivot(reactionsThatProduce.get(j), A);
					if( pivotJ != null )
						k = j;
					else {
						isMinimal = false;
						break;
					}
				}
	 		}
			
			if( isMinimal ) {
				//System.out.println(numRecursions+") Target="+target.getId()+" - #R = "+reactionsThatProduce.size()+". r="+reactionsThatProduce.get(k));
				return reactionsThatProduce.get(k);
			}
		}
		
		return null;
	}
		
	public boolean setHasPrecursors(List<Compound> set)
	{
		if( set == null )
			return false;
		
		for(Compound c: set)
		{
			if( c.isPrecursor() )
				return true;
		}
		return false;
	}

	/* Introducing the next two functions in order to change the merge algorithm in order to choose between
	 * a version that does the minimality test and other version that just checks if the reaction has been analyzed
	 */
	public Reaction findReactionThatProduces(Compound target, HashMap<String, Compound> Hk, boolean allowReactionWithOnlySourcesHkSubstrates)
	{
		List<Reaction> reactions = new ArrayList<Reaction>(target.getReactionsThatProduce(false));
		if( reactions.size() == 0)
			return null;
		
		// prefer reactions that have at least one substrate which is not source nor is in Hk
		for(Reaction r: reactions) {
			Compound pivot = choosePivot(r, Hk);
			if( pivot != null )
				return r;
		}
		return null;
	}
	
	public List<PrecursorSet> precursorSetCartesianUnion(List<PrecursorSet> l1, List<PrecursorSet> l2)
	{
		List<PrecursorSet> c = new ArrayList<PrecursorSet>();
		for(PrecursorSet s1: l1)
		{
			for(PrecursorSet s2: l2)
			{
				PrecursorSet s = new PrecursorSet(s1);
				s.add(s2);
				c.add(s);
			}
		}
		return c;
	}

	
	
	
	static public List<PrecursorSet> reduceToMinimalPrecursorSets(List<PrecursorSet> l)
	{
		// transform all multisets in sets (eliminate duplicated compounds)
		for(int i = 0; i < l.size(); i++)
			transformMultisetInPrecursorSet(l.get(i));
		
		// goes throw the list of solutions and eliminates the solutions that contain other solutions
		List<PrecursorSet> reducedL = new ArrayList<PrecursorSet>(l);
		for(int i = 0; i < l.size(); i++)
		{
			for(int j = 0; j < l.size(); j++)
			{
				if(i == j)
					continue;
				if( !l.get(j).isFlag() && solutionIsSubsetOf(l.get(j), l.get(i), false ) )
				{
					reducedL.remove(l.get(i));
					l.get(i).setFlag( true );
					break;
				}
			}
		}

		// return the resulting set		
		return reducedL;
	}

	/**
	 * Checks if a PrecursorSet s1 is subset of a PrecursorSet s2.
	 * Considers both the list of precursor compounds and the list of Reactions.
	 * 
	 * @param s1
	 * @param s2
	 * @param properSubset
	 * @return true if s2 is a subset of s1
	 */
	static public boolean precursorSetIsSubsetOf(PrecursorSet s1, PrecursorSet s2, boolean properSubset)
	{
		// IMPORTANT: This change has the intention to preserve other solutions different then the empty set.
		/*if( specialEmptySet )
		{
			if( s1.isEmpty() )
				return false;
		}*/
		// MWa delete the following two lines after implementation of own minimization method

		if( s1.isEmpty() )
			return false;

		HashSet<Compound> hashPrecursors = new HashSet<Compound>();
		for(int i = 0; i < s2.getPrecursors().size(); i++)
			hashPrecursors.add(s2.getPrecursors().get(i));
		for(int j = 0; j < s1.getPrecursors().size(); j++)
		{
			if( ! hashPrecursors.contains(s1.getPrecursors().get(j)) )
				return false;
		}
		hashPrecursors.clear();

		if(s1.getPrecursors().size() != s2.getPrecursors().size()){
			return true;
		}

		HashSet<Reaction> hashReactions = new HashSet<Reaction>();
		Reaction[] tmpArray = s2.getReactions().toArray(new Reaction[0]);
		for(int i = 0; i < s2.getReactions().size(); i++)
			hashReactions.add(tmpArray[i]);

		tmpArray = s1.getReactions().toArray(new Reaction[0]);
		for(int j = 0; j < s1.getReactions().size(); j++)
		{
			if( ! hashReactions.contains(tmpArray[j]) )
				return false;
		}

		/*if((s2.getReactions().size() > s1.getReactions().size() )){
			System.out.println(s1.toString() + "\n" + s2.toString());
			System.exit(0);
		}*/

		// check if they are the same...
		return !properSubset || (s2.getReactions().size() > s1.getReactions().size() );  // Works only if are sets (and no multisets)

	}
	
	static public void transformMultisetInPrecursorSet(PrecursorSet s)
	{
		HashSet<Compound> hash = new HashSet<Compound>();
		for(int i = s.getPrecursors().size()-1; i >= 0; i--)
		{
			if( hash.contains(s.getPrecursors().get(i)) )
				s.getPrecursors().remove(i);
			else
				hash.add(s.getPrecursors().get(i));
		}
	}
	
	
	static public boolean solutionIsSubsetOf(PrecursorSet s1, PrecursorSet s2, boolean properSubset)
	{
		// IMPORTANT: This change has the intention to preserve other solutions different then the empty set.
		/*if( specialEmptySet )
		{
			if( s1.isEmpty() )
				return false;
		}*/
		
		HashSet<Compound> hash = new HashSet<Compound>();
		for(int i = 0; i < s2.getPrecursors().size(); i++)
			hash.add(s2.getPrecursors().get(i));
		for(int j = 0; j < s1.getPrecursors().size(); j++)
		{
			if( ! hash.contains(s1.getPrecursors().get(j)) )
				return false;
		}

		// check if they are the same...
		return !properSubset || (s2.getPrecursors().size() > s1.getPrecursors().size());  // Works only if are sets (and no multisets)
	}	
	
	
	/*
	 * The methods below are transforming a "empty set" producer in a set containing only 1 special and new compound to replace the empty set, avoiding
	 * that, by minimalization of the solutions, this "partial" solution is lost.
	 * These methods should work only for the "-mP" method, i.e, by the traversal of the metabolic network, but not yet for the "-mM" merge option.
	 */
	public boolean containsEmptySet(List<PrecursorSet> solutions) 
	{
		for(int i = solutions.size()-1; i >= 0; i--)
		{
			if( solutions.get(i).isEmpty() )
				return true;
		}
		return false;
	}
		
	public void transformEmptySetForTarget(List<PrecursorSet> solutions, Compound c)
	{
		for(int i = solutions.size()-1; i >= 0; i--)
		{
			if( solutions.get(i).isEmpty() )
			{
				Compound emptyProducerForC = createEmptyProducer(c);
				solutions.get(i).addPrecursor(emptyProducerForC);
				break;
			}
		}
	}
	
	public Compound createEmptyProducer(Compound target)
	{
		String precursorKey = target.getId()+" {by empty}";
		if( getNetwork().getCompounds().containsKey(precursorKey) )
			return getNetwork().getCompounds().get(precursorKey);

		// Otherwise, create it...
		Compound emptyCompound = getNetwork().addCompound(precursorKey, precursorKey, "EMPTY SET PRODUCERS");
		emptyCompound.setTopologicalPrecursor(true);
		emptyCompound.setEmptyCompound(true);
		emptyCompound.setCycledCompound(target);
//		Reaction emptyProducer = analyst.getNetwork().addNewReaction("emptyProducerFor"+target.id, null, "emptyProducerFor"+target.id, false);
//		emptyProducer.addProduct(target);
//		emptyProducer.addSubstrate(precursor);
		return emptyCompound;
	}
	
	public void printSolutions(List<PrecursorSet> setOfSets, List<Compound> tcs)
	{
		if( setOfSets == null )
		{
			System.out.println("no solution.");
			return;
		}

		// first, order each solution
		for (PrecursorSet set: setOfSets) 
		{
			Collections.sort(set.getPrecursors(), null);
		}
		
		// then order the set of solutions
		Collections.sort(setOfSets, null);
		
		int countingSets = 1;
		for (PrecursorSet set: setOfSets) 
		{
			String separator = "";
			System.out.println(countingSets++ + ")");
			System.out.print("Precursors : { ");
			for(Compound compound: set.getPrecursors())
			{	
				String id = compound.getId();
				id = StringUtils.sbmlDecode(compound.getId());

				System.out.print(separator + id);
				separator = ", ";
			}
			System.out.println(" }");
			
			System.out.print("Supplementary bootstrap compounds : { ");
			separator = "";
			for(Compound compound: set.getBootstraps())
			{	
				String id = compound.getId();
				id = StringUtils.sbmlDecode(compound.getId());
					
				System.out.print(separator + id);
				separator = ", ";
				
			}
			System.out.println(" }");
			
			if( set.getPrecursors().isEmpty() && set.getReactions() != null ) {
				System.out.println("\n-----\n\nReaction in the cycle:\n");
				for(Reaction r: set.getReactions()) { 
					System.out.println(r.getId());
				}
				System.out.println("\n\n----\n");
			}
			
		}
	}
	
}
