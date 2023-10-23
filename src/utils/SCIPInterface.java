package utils;

import java.util.List;

import metabolicNetwork.Compound;
import application.PrecursorSet;
import de.zib.jscip.nativ.NativeScipException;
import de.zib.jscip.nativ.jni.JniScip;
import de.zib.jscip.nativ.jni.JniScipCons;
import de.zib.jscip.nativ.jni.JniScipConsKnapsack;
import de.zib.jscip.nativ.jni.JniScipConsLinear;
import de.zib.jscip.nativ.jni.JniScipVar;

public class SCIPInterface implements OptimisationInterface {

	@Override
	public void startup() {
		// TODO Auto-generated method stub

	}

	@Override
	public PrecursorSet findNextMinimalPrecursor(List<Compound> sources,
			Compound target, PrecursorSet lastSolution, double bigM,
			double epsilon1) {
		try
		{

			/**@note C pointer are longs in the JNI interface */
			long scip;
			long consLinear;
			long consKnapsack;
			long consLogicor;
			long consSetpart;
			long consSetpack;
			long consSetcover;

			/* create the SCIP environment */
			JniScip env = new JniScip();

			/* create the SCIP variable environment */
			JniScipVar envVar = new JniScipVar();

			/* create the SCIP  constraint environment */
			JniScipCons envCons = new JniScipCons();

			/* create the SCIP knapsack constraint environment */
			JniScipConsKnapsack envConsKnapsack = new JniScipConsKnapsack();

			/* create the SCIP linear constraint environment */
			JniScipConsLinear envConsLinear = new JniScipConsLinear();

			/* create SCIP instance */
			scip = env.create();
		}
		catch (NativeScipException e)
		{
			System.out.println(e.getMessage());
		}
		return null;
	}

	@Override
	public PrecursorSet findNextDuplicatingMachineryMinimalPrecursor(
			List<Compound> sources, Compound target, PrecursorSet lastSolution,
			double bigM, double epsilon1, double epsilon2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public PrecursorSet findNextSteadyStateMinimalPrecursor(
			List<Compound> sources, Compound target, PrecursorSet lastSolution,
			double bigM, double epsilon1) {
		// TODO Auto-generated method stub
		return null;
	}
	
	public List<PrecursorSet> findNextsMinimalPrecursor(List<Compound> sources,
			Compound target, List<PrecursorSet> lastSolutions, double bigM, double epsilon1) {
		return null;
	}
	
	@Override
	public void finish() {
		// TODO Auto-generated method stub

	}

}
