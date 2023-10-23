/* 
PITUFO - A software tool to find all minimal precursor sets for a given set of targets in metabolic networks.

Copyright (C) 2011 Ludovic Cottret (l.cottret@gmail.com <mailto:l.cottret@gmail.com>), Paulo Vieira Milreu  (paulovieira@milreu.com.br <mailto:paulovieira@milreu.com.br>)      
This file is part of PITUFO.

PITUFO is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PITUFO is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PITUFO.  If not, see <http://www.gnu.org/licenses/>.
*/
package metabolicNetwork;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

public class Reaction implements Serializable, Comparable<Reaction> {

	//Auto generated serialVersionUID
	private static final long serialVersionUID = -6342240390701293724L;
	String id;
	String name;
	
	private boolean reversible = false;
	private Reaction reverseDirection = null;
	private int idxBitSet;
	private boolean flag = false;
	
	HashMap<String, Compound> substrates = new HashMap<String, Compound>();
	HashMap<String, Compound> produces = new HashMap<String, Compound>();
	HashMap<String, Compound> sideCompounds = new HashMap<String, Compound>();
	HashMap<String, Compound> cofactors = new HashMap<String, Compound>();
	

	// stoichiometry
	HashMap<String, Double> substratesStoich = new HashMap<String, Double>();
	HashMap<String, Double> productsStoich = new HashMap<String, Double>();
	
	HashMap<String, Reaction> compressedReactions = new HashMap<String, Reaction>();

	public Reaction(String id, String name, boolean reversible) 
	{
		this.id = id;
		this.name = name;
		this.reversible = reversible;		
	}
	
	public void addSubstrate(Compound c)
	{
		if( !substrates.containsKey(c.getId()) )
		{
			substrates.put(c.getId(), c);
			c.addSubstrateOf(this);
		}
	}
	
	public void addSubstrate(Compound c, double stoichiometry){
		addSubstrate(c);
		if( ! substratesStoich.containsKey(c.getId())){
			substratesStoich.put(c.getId(), stoichiometry);
		}
		
	}
	
	public void removeSubstrate(Compound c) {
		c.removeSubstrateOf(this);
		substrates.remove(c.getId());
		substratesStoich.remove(c.getId());
		
		
		//c.getSubstrateOf().remove(this);
	}
	
	public void addSubstrates(List<Compound> cs)
	{
		for(Compound c: cs) 
			addSubstrate(c);
	}	

	public void addProduct(Compound c)
	{
		if( !produces.containsKey(c.getId()) )
		{
			produces.put(c.getId(), c);
			c.addProducedBy(this);
		}
	}
	
	public void addProduct(Compound c, double stoichiometry){
		addProduct(c);
		if( ! productsStoich.containsKey(c.getId())){
			productsStoich.put(c.getId(), stoichiometry);
		}
	}
	
	public void removeProduct(Compound c) {
		c.removeProducedBy(this);
		produces.remove(c.getId());
		productsStoich.remove(c.getId());
		
//		c.getProducedBy().remove(this);
	}
	
	public void addProducts(Collection<Compound> cp)
	{
		for(Compound c: cp)
			addProduct(c);
	}
	
	public boolean ready(List<Compound> inputList, List<Compound> bootstrapList, List<Compound> synthetized)
	{
	  int availableAsSeed = 0;
	  int unavailable = 0;
	  for (Iterator<Compound> iter = substrates.values().iterator(); iter.hasNext();) 
	  {
		  Compound c = (Compound) iter.next();
		  if( inputList.contains(c) || synthetized.contains(c) )
			  availableAsSeed++;
		  else
			  unavailable++;
	  }

	  // If the substrates are OK, the reaction is ready to fire...
	  if( (unavailable == 0) && (availableAsSeed > 0) )
	  	return true;
	  
	  return false;
	}
	
	public boolean canSynthetize(Compound compound)
	{
		return produces.containsKey(compound.getId());
	}
	
	public String toString() {
		StringBuffer str = new StringBuffer(id+": ");
		
		int n = 0; 
		for(Compound cpd: substrates.values()) {
			if(n!=0)  {
				str.append(" + ");
			}
			
			if(substratesStoich.containsKey(cpd.getId())){
				str.append(substratesStoich.get(cpd.getId()) + " ");
			}
			else{
				str.append("1.0 ");
			}
			str.append(cpd.id);
			n++;
		}
		
		if(reversible == true) {
			str.append(" <-> ");
		}
		else {
			str.append(" -> ");
		}
		
		n = 0;
		for(Compound cpd: produces.values()) {
			if(n!=0)  {
				str.append(" + ");
			}
			if(productsStoich.containsKey(cpd.getId())){
				str.append(productsStoich.get(cpd.getId()) + " ");
			}
			else{
				str.append("1.0 ");
			}
			str.append(cpd.id);
			n++;
		}
		return str.toString();
	}
	
	public void setReverse(Reaction r)
	{
		// If there is already some reaction set as reverse, change it to null.
		if( reverseDirection != null )
			reverseDirection.reverseDirection = null;
		
		// Set the new reverse reaction
		reverseDirection = r;
		
		// If there is a reverse reaction defined
		if( reverseDirection != null )
			reverseDirection.reverseDirection = this;
	}
	
	public Reaction getReverseReaction()
	{
		if( reversible )
			return reverseDirection;
		return null;
	}
	
	public boolean isReversible() 
	{
		return reversible || reverseDirection != null;
	}

	public void setReversible(boolean reversible) {
		this.reversible = reversible;
		if( !reversible )
			setReverse(null);
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public HashMap<String, Compound> getSubstrates() {
		return substrates;
	}

	public void setSubstrates(HashMap<String, Compound> substrates) {
		this.substrates = substrates;
	}
	
	public HashMap<String, Double> getSubstratesStoich(){
		return this.substratesStoich;
	}

	public HashMap<String, Compound> getProduces() {
		return produces;
	}

	public void setProduces(HashMap<String, Compound> produces) {
		this.produces = produces;
	}
	
	public HashMap<String, Double> getProductsStoich(){
		return this.productsStoich;
	}

	public HashMap<String, Compound> getSideCompounds() {
		return sideCompounds;
	}

	public void setSideCompounds(HashMap<String, Compound> sideCompounds) {
		this.sideCompounds = sideCompounds;
	}

	public HashMap<String, Reaction> getCompressedReactions() {
		return compressedReactions;
	}

	public void setCompressedReactions(HashMap<String, Reaction> compressedReactions) {
		this.compressedReactions = compressedReactions;
	}
	
	public void initSideCompounds() {
		getSideCompounds().clear();
		for(Compound c: produces.values()) {
			getSideCompounds().put(c.getId(), c);
		}
	}

	public boolean isFlag() {
		return flag;
	}

	public void setFlag(boolean flag) {
		this.flag = flag;
	}

	public void clearFlag() {
		flag = false;
	}
	
	public void setFlag() {
		flag = true;
	}
	
	@Override
	public boolean equals(Object obj) {
		if( !(obj instanceof Reaction) ) {
			return false;
		}
		Reaction other = (Reaction)obj;
		
		if(! id.equals(other.getId())){
			return false;
		}
		
		if( isReversible() != other.isReversible() )
			return false;
		
		// check if the two objects have the same substrates
		if( getSubstrates().size() != other.getSubstrates().size() || getProduces().size() != other.getProduces().size()) {
			return false;
		}
		List<Compound> myList = new ArrayList<Compound>(getSubstrates().values());
		for(Compound c: other.getSubstrates().values()) {
			if( myList.contains(c) && this.getSubstrateStochiometricValue(c) == other.getSubstrateStochiometricValue(c)) {
				myList.remove(c);
			}
		}
		if( myList.size() != 0 ) // if they are equal, I removed all element of the other list
			return false;
		
		// check if the two objects have the same products
		myList = new ArrayList<Compound>(getProduces().values());
		for(Compound c: other.getProduces().values()) {
			if( myList.contains(c) && this.getProductStochiometricValue(c) == other.getProductStochiometricValue(c)) {
				myList.remove(c);
			}
		}
		if( myList.size() != 0 ) // if they are equal, I removed all element of the other list
			return false;
		
		return true;
	}

	public boolean haveSameSubstratesAndProducts(Object obj) {
		if( !(obj instanceof Reaction) ) {
			return false;
		}

		Reaction other = (Reaction)obj;
		boolean sameSubstrates = false;

		// if they do not have the same number of substrates or products
		boolean hasSameNbSubstratesProducts = false;

		if(this.getSubstrates().values().size() == other.getSubstrates().values().size() &&
				this.getProduces().values().size() == other.getProduces().values().size()){
			hasSameNbSubstratesProducts = true;
		}
		if(this.getSubstrates().values().size() == other.getProduces().values().size() &&
				this.getProduces().values().size() == other.getSubstrates().values().size()){
			hasSameNbSubstratesProducts = true;
		}
		if(hasSameNbSubstratesProducts == false){
			return false;
		}

		List<Compound> myList = new ArrayList<Compound>(getSubstrates().values());
		for(Compound c: other.getSubstrates().values()) {
			if( myList.contains(c) && this.getSubstrateStochiometricValue(c) == other.getSubstrateStochiometricValue(c)) {
				myList.remove(c);
			}
		}

		if( myList.size() != 0 ){ // if they are equal, I removed all element of the other list
			// check if the products of 'other' are identical with substrates of 'this'
			myList = new ArrayList<Compound>(getSubstrates().values());
			for(Compound c: other.getProduces().values()) {
				if( myList.contains(c) && this.getSubstrateStochiometricValue(c) == other.getProductStochiometricValue(c)) {
					myList.remove(c);
				}
			}

			if( myList.size() != 0 ){// neither substrates nor products of 'other' are identical with substrates of 'this' 
				return false;
			}
		}

		else{
			sameSubstrates = true;
		}

		boolean sameProducts = false;
		// check if the two objects have the same products
		myList = new ArrayList<Compound>(getProduces().values());
		if(sameSubstrates){ // substrates are equal -> compare products
			for(Compound c: other.getProduces().values()) {
				if( myList.contains(c) && this.getProductStochiometricValue(c) == other.getProductStochiometricValue(c)) {
					myList.remove(c);
				}
			}

			if(myList.size() == 0){
				sameProducts = true;
			}
		}
		else{ // compare products of 'this' with substrates of 'other'
			for(Compound c: other.getSubstrates().values()) {
				if( myList.contains(c) && this.getProductStochiometricValue(c) == other.getSubstrateStochiometricValue(c)) {
					myList.remove(c);
				}
			}
		}

		if( myList.size() != 0 ) // if they are equal, I removed all element of the other list
			return false;

		if(this.isReversible() || (sameSubstrates && sameProducts && ! other.isReversible())){
			return true;
		}
		return false;
	}

	public double getProductStochiometricValue(Compound c) {
		if (productsStoich.containsKey(c.getId()))
			return productsStoich.get(c.getId());
		return 0;
	}	
	
	public double getSubstrateStochiometricValue(Compound c) {
		if (substratesStoich.containsKey(c.getId()))
			return substratesStoich.get(c.getId());
		return 0;
	}	
	
	public Collection<Compound> getCofactors() {
		return this.cofactors.values();
	}

	public void setCofactors(HashMap<String, Compound> cofactors) {
		this.cofactors = cofactors;
	}
	
	public void addCofactor(Compound c){
		this.cofactors.put(c.getId(), c);
	}
	
	public void addCofactors(Collection<Compound> collection){
		for(Compound c: collection){
			this.cofactors.put(c.getId(), c);
		}
	}

	public int getIdxBitSet() {
		return idxBitSet;
	}

	public void setIdxBitSet(int idxBitSet) {
		this.idxBitSet = idxBitSet;
	}

	@Override
	public int compareTo(Reaction o) {
		return this.toString().compareTo(o.toString());
	}

}
