package utils;

import ilog.concert.IloAddable;
import ilog.concert.IloException;
import ilog.concert.IloIntVar;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloRange;
import ilog.cplex.IloCplex;
import ilog.cplex.IloCplexModeler;


public class SasitaCPLEXModelHolder {
	public static enum SasitaModelType {NORMAL, DUPMACH,STEADYSTATE};

	private IloCplex cplex;
	private IloNumVar[] varX = null;
	private IloIntVar[] varInd = null;
	private IloCplexModeler modeler = new IloCplexModeler();
	private SasitaModelType modelType = null;
	private boolean isModelSet = false;
	private int sizeOfTheBiggestSolutionFound = 0;
	private IloRange sizeConstraint = null;
	private IloLinearNumExpr fobj = null;
	
	public boolean isModelSet() {
		return isModelSet;
	}

	public IloCplex getCplex() {
		return cplex;
	}
	
	public IloNumVar[] getVarX() {
		return varX;
	}
	
	public IloIntVar[] getVarInd() {
		return varInd;
	}
	
	public SasitaModelType getModelType() {
		return modelType;
	}
	
	public int getSizeOfBiggestSolution() {
		return sizeOfTheBiggestSolutionFound;
	}

	public void addSolutionSizeConstraint(int lastSizeConstraint) {
		if(lastSizeConstraint > this.sizeOfTheBiggestSolutionFound){
			this.sizeOfTheBiggestSolutionFound = lastSizeConstraint;
			try{
				if (this.sizeConstraint != null){
					cplex.remove(this.sizeConstraint);
				}
				this.sizeConstraint = modeler.addRange(this.sizeOfTheBiggestSolutionFound,this.fobj,this.varInd.length,"FObjInfBound_>=" + this.sizeOfTheBiggestSolutionFound);
				this.cplex.add(this.sizeConstraint);
			} catch (IloException e) {
				e.printStackTrace();
				System.err.println("Error in the CPLEX Model Holder, Aborting.");
				System.exit(-1);
			}
		}
	}

	public void setModel(IloCplex cplex, IloNumVar[] varX, IloIntVar[] varInd, SasitaModelType modelType){
		this.cplex = cplex;
		this.varX = varX;
		this.varInd = varInd;
		this.modelType = modelType;
		this.isModelSet = true;
		try {
			this.fobj = this.modeler.linearNumExpr();
			for (IloIntVar Ii : varInd) {
				this.fobj.addTerm(1.0, Ii);
			}
		} catch (IloException e) {
			e.printStackTrace();
			System.err.println("Error in the CPLEX Model Holder, Aborting.");
			System.exit(-1);
		}
	}
	
	/*
	 * Sets a start point where every variable is 1.0,
	 * for a fresh resolve.
	 */
	public void touchModelToResolve(){
		try {
			int nr = this.varX.length;
			int nind = this.varInd.length;
			double[] rv = new double[nr];
			double[] indv = new double[nind];
			for(int i = 0; i<nr;++i){
				rv[i] = 1.0;
			}
			for(int i = 0; i<nind;++i){
				indv[i] = 1.0;
			}
			this.cplex.setStart(rv, null, this.varX, null, null, null);
			this.cplex.setStart(indv, null, this.varInd, null, null, null);
			this.cplex.addMIPStart(varInd, indv);
			this.cplex.addMIPStart(varX, rv);
		} catch (IloException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public void clearModel(){
		try {
			this.cplex.clearModel();
		} catch (IloException e) {
			e.printStackTrace();
			System.err.println("Error in the CPLEX Model Holder, Aborting.");
			System.exit(-1);
		}
		this.cplex.end();
		this.cplex = null;
		this.varX = null;
		this.varInd = null;
		this.modelType = null;
		this.fobj = null;
		this.sizeConstraint = null;
		this.sizeOfTheBiggestSolutionFound = 0;
		this.isModelSet = false;
	}
}
