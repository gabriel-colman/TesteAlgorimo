package utils;

import java.util.List;

import application.PrecursorSet;
import metabolicNetwork.Compound;

public interface OptimisationInterface {
	void startup();

	PrecursorSet findNextMinimalPrecursor(List<Compound> sources,
			Compound target, PrecursorSet lastSolution, double bigM, double epsilon1);
	
	PrecursorSet findNextDuplicatingMachineryMinimalPrecursor(
			List<Compound> sources, Compound target, PrecursorSet lastSolution, double bigM, double epsilon1, double epsilon2);
	
	PrecursorSet findNextSteadyStateMinimalPrecursor(
			List<Compound> sources, Compound target, PrecursorSet lastSolution, double bigM, double epsilon1);

	List<PrecursorSet> findNextsMinimalPrecursor(List<Compound> sources,
			Compound target, List<PrecursorSet> lastSolutions, double bigM, double epsilon1);
		
	void finish(); 
}
