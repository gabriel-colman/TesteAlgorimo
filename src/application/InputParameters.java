package application;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import metabolicNetwork.Compound;

public class InputParameters {
	private static List<Compound> inputCompounds = new ArrayList<Compound>();
	private static List<Compound> bootstrapCompounds = new ArrayList<Compound>();
	private static List<Compound> userDefinedPrecursors = new ArrayList<Compound>();
	private static List<Compound> targetCompounds = new ArrayList<Compound>();
	public static List<Compound> splitSources = new ArrayList<Compound>();

	public static boolean decodeSbml = true;
	public static boolean printForward = false;
	public static boolean printNetwork = false;
	public static boolean oneByOne = true;
	public static boolean keepForwardResultInTheTree = false;
	public static boolean minimalityCheck = true;
	public static boolean mergeOnlyWithMinimalReactions = true;
	public static boolean forwardPropagation = false;
	public static boolean reductionThroughMaximalTarget = false;
	public static boolean precursorIfProducedOnlyByReversible = true;
	public static boolean checkSolutions = false;
	public static boolean mergeKeepingSideCompounds = true;
	public static int     maxTimePerTarget = 0;//1 * 60 * 1000;
	public static boolean emptyCyclesEnumeration = false;
	public static boolean emptyCyclesCut = false;
	public static boolean randomChoices = true;
	public static int     stopNoNew = 0;
	public static int	  nbThreads = 1;
	public static String  forceSources = "";
	public static String  filterPairedCofactors = "";
	public static File  dirPrecursorSolutions = new File("");
	public static boolean considerOnlyUserDefinedPrecursors = false;
	public static boolean stopInByProducts = false;
	public static boolean verbose = false;
	public static boolean noPreprocessing = false;
	public static double epsilon1 = 0.1;
	public static double epsilon2 = 0.1;
	public static double bigM = 1000;
	public static boolean modeDuplicatingMachinery = false;
	public static boolean modeSteadyState = false;
	public static boolean modeCheckSolutionsOnly = false;
	public static boolean checkMinimality;
	public static boolean addReactionsToSolutions = false;
	public static boolean addCumulatedCompoundsToSolutions = false;
	public static boolean fbaCheck;
	public static int TiLim;

	
	public static List<Compound> getInputCompounds() {
		return inputCompounds;
	}
	public static void setInputCompounds(List<Compound> inputCompounds) {
		InputParameters.inputCompounds = inputCompounds;
	}
	public static List<Compound> getUserDefinedPrecursors() {
		return userDefinedPrecursors;
	}
	public static void setUserDefinedPrecursors(List<Compound> userDefinedPrecursors) {
		InputParameters.userDefinedPrecursors = userDefinedPrecursors;
	}
	public static List<Compound> getBootstrapCompounds() {
		return bootstrapCompounds;
	}
	public static void setBootstrapCompounds(List<Compound> bootstrapCompounds) {
		InputParameters.bootstrapCompounds = bootstrapCompounds;
	}
	public static List<Compound> getTargetCompounds() {
		return targetCompounds;
	}
	public static void setTargetCompounds(List<Compound> targetCompounds) {
		InputParameters.targetCompounds = targetCompounds;
	}
	public static boolean getConsiderOnlyUserDefinedPrecursors() {
		return considerOnlyUserDefinedPrecursors;
	}
	public static void setConsiderOnlyUserDefinedPrecursors(
			boolean considerOnlyUserDefinedPrecursors) {
		InputParameters.considerOnlyUserDefinedPrecursors = considerOnlyUserDefinedPrecursors;
	}

}
