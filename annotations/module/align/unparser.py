#
# Contain unparsing procedures.
#

import sys

#-------------------------------------------

def unparseToC(vars, annot_body_code, indent, extra_indent):
    '''Unparse to C/C++ code'''

    if len(vars) == 0:
        return annot_body_code

    s = '\n'
    s += indent + '#pragma disjoint ('
    for i, v in enumerate(vars):
        if i > 0:
            s += ', '
        s += '*' + __printAddressC(v.var_name, v.dimensions)
    s += ') \n'
    s += indent + 'if ((('
    for i, v in enumerate(vars):
        if i > 0:
            s += '|'
        s += '(int)(' + __printAddressC(v.var_name, v.dimensions) + ')'
    s += ') & 0xF) == 0) {\n'
    for v in vars:
        s += indent + extra_indent 
        s += '__alignx(16,' + __printAddressC(v.var_name, v.dimensions) + ');\n'
    s += annot_body_code.replace('\n', '\n' + extra_indent)
    s += '\n'
    s += indent + '} else {\n'
    s += annot_body_code.replace('\n', '\n' + extra_indent)
    s += '\n'
    s += indent + '}\n'
    s += indent
    return s

#-------------------------------------------

def unparseToFortran(vars, annot_body_code, indent, extra_indent):
    '''Unparse to Fortran code'''

    print 'error: Fortran is not yet supported in alignment module'
    sys.exit(1)

#-------------------------------------------

def __printAddressC(var_name, dimensions):
    '''Return the starting address location of the given variable (in C/C++)'''

    dimensions = dimensions[:]
    dimensions.remove(None)
    s = str(var_name)
    if len(dimensions) > 0:
        s += '['
        s += ']['.join(map(str, dimensions))
        s += ']'
    return s

#-------------------------------------------

def __printAddressFortran(var_name, dimensions):
    '''Return the starting address location of the given variable (in Fortran)'''

    dimensions = dimensions[:]
    dimensions.remove(None)
    s = str(var_name)
    if len(dimensions) > 0:
        s += '('
        s += ','.join(map(str, dimensions))
        s += ')'
    return s        
    
