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
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

public class Compound implements Comparable<Compound>, Serializable {
	
	
	//Auto generated serialVersionUID
	private static final long serialVersionUID = -6326653177519060454L;
	String id;
	String name;
	String compartment;
	private HashMap<String, Integer> atoms;
	private boolean isBoundary = false;
	private boolean userDefinedPrecursor = false;
	private boolean topologicalPrecursor = false;
	private boolean allowed = true;
	private boolean bootstrap = false;
	private boolean target = false;
	private boolean emptyCompound = false;
	private boolean flag = false;
	private boolean forcedSource = false;
	private int idxBitSet;
	private boolean cycleFlag = false;
	private boolean positiveByProduct = false;
	private boolean highDegree = false;
	
	private Compound cycledCompound;
	private boolean markAsAlreadyVisitedAtSameLevel = false;
	
	List<Reaction> producedBy = new ArrayList<Reaction>();
	List<Reaction> substrateOf = new ArrayList<Reaction>();
	

	// for tarjan stronly connected algorithm
	private Integer tarjanIndex = -1;
	private Integer tarjanLowlink = -1;
	
	public Compound(String id, String name, String compartment) {
		this.id = id;
		this.name = name;
		this.compartment = compartment;
	}
	
	public String getCompartment() {
		return compartment;
	}

	public void setCompartment(String compartment) {
		this.compartment = compartment;
	}

	public void setFlag(boolean flag) {
		this.flag = flag;
	}

	
	public void copyPropertiesFrom(Compound source) {
		this.userDefinedPrecursor = source.userDefinedPrecursor;
		this.topologicalPrecursor = source.topologicalPrecursor;
		this.allowed = source.allowed;
		this.bootstrap = source.bootstrap;
		this.target = source.target;
		this.flag = source.flag;
		this.emptyCompound = source.emptyCompound;
		this.idxBitSet = source.idxBitSet;
		this.highDegree = source.highDegree;
		this.isBoundary = source.isBoundary();
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

	public List<Reaction> getProducedBy() {
		return producedBy;
	}

	public void setProducedBy(List<Reaction> producedBy) {
		this.producedBy = producedBy;
	}

	public List<Reaction> getSubstrateOf() {
		return substrateOf;
	}

	public void setSubstrateOf(List<Reaction> substrateOf) {
		this.substrateOf = substrateOf;
	}

	public Compound getCycledCompound() {
		return cycledCompound;
	}

	public void setCycledCompound(Compound cycledCompound) {
		this.cycledCompound = cycledCompound;
	}

	

	@Override
	public String toString() {
		
		return "<" + id +">";
	}
		
	public int compareTo(Compound c) {
		return toString().compareTo(c.toString());
	}
	
	public void addProducedBy(Reaction r)
	{
		if( !producedBy.contains(r))
		{
			producedBy.add(r);
		}
	}
	
	public void removeProducedBy(Reaction r){
		producedBy.remove(r);
	}
	
	public void removeSubstrateOf(Reaction r){
		substrateOf.remove(r);
		
	}
	
	public void addSubstrateOf(Reaction r)
	{
		if( !substrateOf.contains(r))
		{
			substrateOf.add(r);
		}
	}

	public boolean isPrecursor() {
		return userDefinedPrecursor || topologicalPrecursor;
	}

	public List<Reaction> getReactionsThatProduce(boolean considerReversibility)
	{
		if( !considerReversibility )
			return new ArrayList<Reaction>(producedBy);
		else
		{
			List<Reaction> rs = new ArrayList<Reaction>(producedBy);
			for(int i = 0; i < substrateOf.size(); i++)
			{
				Reaction r = substrateOf.get(i);
				if( r.isReversible() )
				   rs.add(r);
			}
			return rs;
		}
				
	}

	public List<Reaction> getReactionsThatConsume(boolean considerReversibility)
	{
		if( !considerReversibility )
			return new ArrayList<Reaction>(substrateOf);
		else
		{
			List<Reaction> rs = new ArrayList<Reaction>(substrateOf);
			for(int i = 0; i < producedBy.size(); i++)
			{
				Reaction r = producedBy.get(i);
				if( r.isReversible() )
				   rs.add(r);
			}
			return rs;
		}
				
	}

	public void setUserDefinedPrecursor(boolean userDefinedPrecursor) 
	{
		if( !isTarget() )
			this.userDefinedPrecursor = userDefinedPrecursor;
	}

	public boolean isUserDefinedPrecursor() {
		return userDefinedPrecursor;
	}

	public boolean isTopologicalPrecursor() {
		return topologicalPrecursor;
	}

	public void setTopologicalPrecursor(boolean topologicalPrecursor) 
	{
		if( !isTarget() )
			this.topologicalPrecursor = topologicalPrecursor;
	}

	public boolean isBoundary() {
		return isBoundary;
	}

	public void setBoundary(boolean isBoundary) {
		this.isBoundary = isBoundary;
	}

	public void setAllowed(boolean allowed) {
		this.allowed = allowed;
	}

	public boolean isAllowed() {
		return allowed;
	}

	public void setBootstrap(boolean bootstrap) {
		this.bootstrap = bootstrap;
	}

	public boolean isBootstrap() {
		return bootstrap;
	}

	public void setTarget(boolean target) 
	{
		this.target = target;
	}

	public boolean isTarget() {
		return target;
	}
	
	public boolean isEmptyCompound() {
		return emptyCompound;
	}

	public void setEmptyCompound(boolean emptyCompound) {
		this.emptyCompound = emptyCompound;
	}

	public void clearFlag() {
		flag = false;
	}
	
	public boolean isFlag() {
		return flag;
	}
	
	@Override
	public boolean equals(Object other) {
		if (other instanceof Compound){
			if ( ((this.id==null && ((Compound)other).id==null) || id.equals(((Compound)other).id))
					&& ((this.name==null && ((Compound)other).name==null) || name.equals(((Compound)other).name))){
				return true;
			}
		}
		return  false;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		if (this.id != null){
			result = prime * result + id.hashCode();
		}
		if (this.name != null){
			result = prime * result + name.hashCode();
		}
		return result;
	}

	public int getTarjanIndex() {
		return tarjanIndex;
	}

	public void setTarjanIndex(int tarjanIndex) {
		this.tarjanIndex = tarjanIndex;
	}

	public int getTarjanLowlink() {
		return tarjanLowlink;
	}

	public void setTarjanLowlink(int tarjanLowlink) {
		this.tarjanLowlink = tarjanLowlink;
	}

	public boolean isCycleFlag() {
		return cycleFlag;
	}

	public void setCycleFlag(boolean cycleFlag) {
		this.cycleFlag = cycleFlag;
	}

	public boolean isMarkAsAlreadyVisitedAtSameLevel() {
		return markAsAlreadyVisitedAtSameLevel;
	}

	public void setMarkAsAlreadyVisitedAtSameLevel(
			boolean markAsAlreadyVisitedAtSameLevel) {
		this.markAsAlreadyVisitedAtSameLevel = markAsAlreadyVisitedAtSameLevel;
	}

	public boolean isForcedSource() {
		return forcedSource;
	}

	public void setForcedSource(boolean forcedSource) {
		this.forcedSource = forcedSource;
	}

	public boolean isPositiveByProduct() {
		return positiveByProduct;
	}

	public void setPositiveByProduct(boolean positiveByProduct) {
		this.positiveByProduct = positiveByProduct;
	}

	public HashMap<String, Integer> getAtoms() {
		return atoms;
	}
	
	public int getNumberOfAtom(String atom){
		if(this.atoms != null && this.atoms.containsKey(atom)){
			return this.atoms.get(atom);
		}
		else{
			return 0;
		}
	}

	public void setAtoms(HashMap<String, Integer> atoms) {
		this.atoms = atoms;
	}

	public int getIdxBitSet() {
		return idxBitSet;
	}

	public void setIdxBitSet(int idxBitSet) {
		this.idxBitSet = idxBitSet;
	}

	public boolean hasHighDegree() {
		return highDegree;
	}

	public void setHighDegree(boolean highDegree) {
		this.highDegree = highDegree;
	}

	
}


class CompoundComparator implements Comparator<Object>{
	@Override
	public int compare(Object o1, Object o2) {
		Compound c1 = (Compound) o1;
		Compound c2 = (Compound) o2;
		
		if(c1.isPrecursor() == c2.isPrecursor()){
			return 0;
		}
		else if(c1.isPrecursor()){
			return 1;
		}
		else{
			return -1;
		}	
	}
}
