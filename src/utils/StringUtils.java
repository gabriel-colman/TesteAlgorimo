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

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class StringUtils {
	
	public static String htmlEncode(String in) {
		String out;
		
		out = in;
		out = out.replaceAll("<[^>]*>", "");
		out = out.replaceAll("&lt;", "less_than");
		out = out.replaceAll("&gt;", "greater_than");
		out = out.replaceAll("&", "&amp;");
		out = out.replaceAll("\"", "&quot;");
		out = out.replaceAll("'", "&apos;");
		
		return out;
	}
	
	public static String sbmlEncode(String in) {
		String out;
		
		out = htmlEncode(in);
		String REGEX = "^\\d";
		
		Pattern pattern = Pattern.compile(REGEX);
        Matcher matcher = pattern.matcher(out);
        
        if(matcher.find()) {
			out = "_".concat(out);
		}
        
        REGEX = "[^0-9A-Za-z_]";
        
        pattern = Pattern.compile(REGEX);
        matcher = pattern.matcher(out);
        
        while(matcher.find()) {
        	String specialCharacter = matcher.group(0);
        	Integer value = specialCharacter.codePointAt(0);
        	String code = "__"+value+"__";
        	
        	out = out.replace(specialCharacter,code);
        }
		
		return out;
	}
	
	public static String sbmlDecode(String in) {
		String out = new String(in);
		
		String REGEX = ".*__(\\d+)__.*";
		
		Pattern pattern = Pattern.compile(REGEX);
        Matcher matcher = pattern.matcher(out);
        
        while(matcher.find()) {
        	String str = matcher.group(1);
        	
          	int codesInt[] = new int[1];
        	codesInt[0] = new Integer(str);
        	
        	String specialCharacter = new String(codesInt,0, codesInt.length);
        	
        	out = out.replace("__"+str+"__", specialCharacter);
        	
        	matcher = pattern.matcher(out);
        }
        
        REGEX = "^_(\\d*).*";
        
        pattern = Pattern.compile(REGEX);
        matcher = pattern.matcher(out);
        
        if(matcher.find()) {
        	String str = matcher.group(1);   	
        	out = out.replaceFirst("^_"+str, str);
        }
		
        return out;
	}
	
	public static void main(String[] args) {
		
//		String in = "<i>cis</i>-zeatin & biosynthesis &lt; 1 &gt; 1; \" ,";
//		
//		System.err.println(htmlEncode(in));
//		
//		String in2 = "2-CYSTEIN rED 4";
//		
//		System.err.println(sbmlEncode(in2));
//		
//		String in = "0.5d0";
//		String in2 = ".5";
//		String in3 = "d1.234f05";
//		
//		System.out.println(transformStoi(in));
//		System.out.println(transformStoi(in2));
//		System.out.println(transformStoi(in3));
		String in = "__124__Pi__124__";
		
		System.out.println(sbmlDecode(in));		
	}
	
	/**
	 * Transforms a stoechiometric coefficient to be compatible with the
	 * SBML annotations
	 */
	public static String transformStoi (String st) {
		
		if(st == null) {
			return "1";
		}
		
		if(st.matches("^\\d$"))
				return st;
		
		String REGEX = "[^\\d]*(\\d*\\.\\d+)[^\\d]*.*";
		
		Pattern pattern = Pattern.compile(REGEX);
        Matcher matcher = pattern.matcher(st);
		
        if(matcher.find()) {
        	String out = matcher.group(1);
        	return out;
        }
        else {
        	return "1";
        }
	}
}