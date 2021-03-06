# Tunguska by Team210 - 64k intro by Team210 at Solskogen 2k19
# Copyright (C) 2018  Alexander Kraus <nr4@z10.info>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import GLSLLexer130

def compress(source):
    #return source
    
    identifier_list = []
    small_identifiers = [ chr(ord('a') + i) for i in range(26) ]
    small_identifiers += [ chr(ord('A') + i) for i in range(26) ]
    small_identifiers += [ a + b for a in small_identifiers for b in small_identifiers ]
    print(small_identifiers)
    
    newcode = ""
    lexer = GLSLLexer130.GLSLLexer130(source)
    token = lexer.token()
    while token != None:
        uniform = False
        
        if (token.tokenName != "SINGLELINE_COMMENT") and (token.tokenName != "MULTILINE_COMMENT") and (token.tokenName != "IDENTIFIER"):
            newcode += token.tokenData
        if token.tokenName in [ "VOID", "FLOAT", "VEC2", "VEC3", "VEC4", "MAT2", "MAT3", "MAT4", "SAMPLER2D", "UNIFORM", "IN_QUALIFIER", "OUT_QUALIFIER", "INOUT_QUALIFIER", "VOID", "VERSION_DIRECTIVE", "DEFINE_DIRECTIVE", "CONST", "INT", "ELSE", "RETURN" ]:
            newcode += " "
        #if token.tokenName == "IDENTIFIER": # No minification
            #newcode += token.tokenData
        if token.tokenName == "IDENTIFIER":
            identifier = token.tokenData
            if not identifier in identifier_list:
                identifier_list += [ identifier ]
            ind = identifier_list.index(identifier) 
            print("ID: ", token.tokenData, ind)
            newcode += small_identifiers[ind]
        token = lexer.token()
    return newcode
