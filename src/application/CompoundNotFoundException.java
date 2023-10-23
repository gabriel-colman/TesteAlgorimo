package application;

import metabolicNetwork.Compound;

public class CompoundNotFoundException extends Exception {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
public 	CompoundNotFoundException(Compound c){
	super(c.getId() + " not found");
}

}
