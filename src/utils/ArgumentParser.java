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
package utils;

/**
 * @author Whit Stockwell, adapted by Ludo COTTRET
 * 
 */
import java.util.*;
public class ArgumentParser {
    public ArgumentParser(String[] args) {
        for (int i = 0; i < args.length; i++) {
            if (args[i].startsWith("-") || args[i].startsWith("/")) {
                int loc = args[i].indexOf("=");
                String key = (loc > 0) ? args[i].substring(1, loc) :
args[i].substring(1);
                String value = (loc > 0) ? args[i].substring(loc+1) :
"";
                options.put(key, value);
            }
            else {
                params.addElement(args[i]);
            }
        }
    }

    public boolean hasOption(String opt) {
        return options.containsKey(opt);
    }

    public String getOption(String opt) {
        return (String) options.get(opt);
    }

    public String nextParam() {
        if (paramIndex < params.size()) {
            return (String) params.elementAt(paramIndex++);
        }
        return null;
    }
    private Vector<String> params = new Vector<String>();
    private Hashtable<String,String> options = new Hashtable<String,String>();
    private int paramIndex = 0;
}

