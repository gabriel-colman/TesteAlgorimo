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
package application;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import pitufo.PrecursorFinder;


/**
 * 
 * This class represents a set of solutions for a given target. 
 *
 */
public class PrecursorSets {
	
	Set<PrecursorSet> solutions = new HashSet<PrecursorSet>();

	boolean changed = false;

	public boolean isChanged() {
		return changed;
	}

	public void setChanged(boolean changed) {
		this.changed = changed;
	}

	/*
	 * addSolution adds a new precursor set to the set of solutions, if it is a new solution (not previously present in the set)
	 * 
	 * it returns true in the case s is a new solution and false otherwise
	 * 
	 */
	public boolean addSolution(PrecursorSet s) {
		for(PrecursorSet sol: solutions) {
			if( sol.equals(s) ) {
				return false;
			}
		}
		setChanged(true);
		solutions.add(s);
		return true;
	}
	
	/*
	 * addSolutions adds all precursor sets in the sols collection, it returns how many new solutions were found.
	 * 
	 */
	public int addSolutions(Collection<PrecursorSet> sols) {
		int numNew = 0;
		for(PrecursorSet sol: sols) {
			if( addSolution(sol) )
				numNew++;
		}
		return numNew;
	}
	
	
	public List<PrecursorSet> minimalize() {
		List<PrecursorSet> sols = new ArrayList<PrecursorSet>(solutions);
		return PrecursorFinder.reduceToMinimalPrecursorSets(sols);
	}
	
	@Override
	public boolean equals(Object obj) {
		if( !(obj instanceof PrecursorSets) ) {
			return false;
		}
		PrecursorSets other = (PrecursorSets)obj;
		if( solutions.size() != other.solutions.size()) {
			return false;
		}
		List<PrecursorSet> myList = new ArrayList<PrecursorSet>(solutions);
		for(PrecursorSet ps: other.solutions) {
			if( myList.contains(ps) ) {
				myList.remove(ps);
			}
		}
		return myList.size() == 0; // if they are equal, I removed all element of the other list
	}
	
}
