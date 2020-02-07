#!/usr/bin/env python

import os
import sys

try:
    from pygccxml import *
except Exception as e:
    print
    print e
    print
    print 'On SL6 try: setup pygccxml v1_9_1  -q p2714b; setup castxml v0_00_00_f20180122'
    print 'On SL7 try: setup pygccxml v1_9_1a -q p2715a; setup castxml v0_00_00_f20180122'
    print
    sys.exit(1)

# Types that ROOT will understand (with a little nudging in some cases)
fundamental_types = ['int', 'float', 'double', 'bool', 'unsigned int',
                     'short', 'short int', 'short unsigned int',
                     'long', 'unsigned long', 'long unsigned int',
                     'long long int', 'unsigned char',
                     'size_t',
                     'std::string']

enums = []

def no_namespace(type):
    return type[type.rfind(':')+1:]


def is_fundamental(type):
    return type in fundamental_types or no_namespace(type) in enums


def translate_type(type):
    # enums are small, right?
    if type in enums: return 'short unsigned int'
    # TTree can't handle long long at all? This will lose information though...
    if type == 'long long int': return 'long'

    return type


def is_vector(type):
    return type[:12] == 'std::vector<'


def vector_contents(type):
    assert is_vector(type)
    return type[12:type.find(',')]


def is_array(type):
    return '[' in type


def array_type(type):
    assert is_array(type)
    # There always seems to be a space before the dimensions
    return type[:type.find('[')-1]


def array_size(type):
    assert is_array(type)
    return int(type[type.find('[')+1:-1])


def type_to_flat_type(type):
    assert not is_fundamental(no_namespace(type))
    assert not is_vector(type)

    if type == 'StandardRecord': return 'FlatRecord'

    return 'Flat'+no_namespace(type)


def base_class(klass):
    assert len(klass.bases) < 2, 'Support for multiple base classes not implemented'
    if len(klass.bases) == 1:
        return klass.bases[0].related_class
    return None


# Recurse to find member variables including all base classes
def variables_inc_bases(klass):
    vars = list(klass.variables())
    # Only accept direct members of the class, not members of any nested class
    vars = [v for v in vars if v.parent == klass]
    base = base_class(klass)
    if base: vars += variables_inc_bases(base)
    return vars


def fq_type(klass):
    # HACK special case for the weights map
    if str(klass.name) == 'Pair':
        return 'std::pair<std::string, std::vector<float>>'

    # HACK we had this within event:: so we could find it, but in this context
    # it's the original ROOT one we want.
#    if str(klass.name) == 'TVector3': return 'TVector3'

    if not klass.parent or str(klass.parent.name) == '::':
        return str(klass.name)
    else:
        return fq_type(klass.parent)+'::'+str(klass.name)


def branch_cmd(name, treeName = 'tr'):
    return 'branch('+treeName+', prefix+"'+name+'", &'+name+', policy);'

def branch_cmd_no_policy(name, treeName = 'tr'):
    return 'branch('+treeName+', prefix+"'+name+'", &'+name+', 0);'



if len(sys.argv) < 1 or len(sys.argv) > 3:
    print 'Usage: gen_flatrecord.py [/path/to/header/outputs/] [/path/to/cxx/outputs/]'
    sys.exit(1)

headerDir = os.getcwd()
if len(sys.argv) >= 2: headerDir = sys.argv[1]
cxxDir = os.getcwd()
if len(sys.argv) >= 3: cxxDir = sys.argv[2]

# Locate the castxml executable
generator_path, generator_name = utils.find_xml_generator()

# Figure out where our source files are
context = os.environ['SBNCODE_DIR']

path = [context, os.environ['ROOT_INC']]

config = parser.xml_generator_configuration_t(
    xml_generator_path=generator_path,
    xml_generator=generator_name,
    include_paths=path,
    cflags='-std=c++1z -DGEN_FLATRECORD_CONTEXT -Wno-unknown-warning-option'#,
#    start_with_declarations='caf::StandardRecord'
    )

print 'Reading from', context+'/sbncode/StandardRecord/StandardRecord.h'
decls = parser.parse([context+'/sbncode/StandardRecord/StandardRecord.h'],
                     config)

global_namespace = declarations.get_global_namespace(decls)
ns = global_namespace.namespace('caf')

enums += [e.name for e in ns.enumerations()]
#enums += ['Experiment'] # isn't currently within the event:: namespace

# Keep track of which classes we've written out so far, for purposes of
# dependency tracking.
emitted = []


disclaimer = '''// This file was auto-generated by gen_flatrecord.py.
// DO NOT EDIT IT DIRECTLY.
// For documentation of the fields see the regular StandardRecord.h'''

# We write individual headers so that the preprocessor will handle our
# dependency order, but put all the code in one .cxx becase it's much faster to
# compile that way.
sys.stdout = file(cxxDir+'/FlatRecord.cxx', 'wa')
print disclaimer
print
print '#include "sbncode/StandardRecord/StandardRecord.h"'
print
print '#include "sbncode/FlatMaker/TreeUtils.h"'
print
print '#include "TTree.h"'

for klass in ns.classes():
    pt = type_to_flat_type(klass.name)

    sys.stdout = file(headerDir+'/'+pt+'.h', 'w')

    class TypeName:
        def __init__(self, t, n, length = -1):
            self.type = t
            self.name = n
            self.length = length

        def is_array(self): return self.length > 0

    basic_mems = []
    class_mems = []
    vec_mems = []
    for v in variables_inc_bases(klass):
        type = str(v.decl_type)

        if is_fundamental(type):
            basic_mems += [TypeName(translate_type(type), v.name)]
        elif is_array(type):
            innert = array_type(type)
            size = array_size(type)
            if is_fundamental(innert):
                vec_mems += [TypeName(innert, v.name, size)]
            else:
                vec_mems += [TypeName(type_to_flat_type(innert), v.name, size)]
        elif is_vector(type):
            innert = vector_contents(type)
            if is_fundamental(innert):
                vec_mems += [TypeName(innert, v.name)]
            else:
                vec_mems += [TypeName(type_to_flat_type(innert), v.name)]
        else:
            if type == 'TVector3':
                sys.stderr.write('Warning: skipping TVector3 '+tn.name+'\n')
                continue
            class_mems += [TypeName(type_to_flat_type(type), v.name)]

    print disclaimer
    print
    print '#pragma once'
    print
#    print '#include "StandardRecord/SREnums.h"'
#    print
    incs = set()
    for tn in class_mems: incs.add(tn.type)
    for tn in vec_mems:
        if not is_fundamental(tn.type):
            incs.add(tn.type)
    for i in incs:
        print '#include "sbncode/FlatMaker/'+i+'.h"'
    print
    print '#include <string>'
    print
    print 'class TTree;'
    print 'namespace caf{class '+klass.name+';}'
    print
    print 'namespace flat'
    print '{'
    print 'class IBranchPolicy;'
    print
    print '/// Flat encoding of \\ref', klass.name
    print 'class', pt
    print '{'
    print 'public:'
    print '  '+pt+'(const std::string& prefix, TTree* tr, const IBranchPolicy* policy);'
    print '  ~'+pt+'();'
    print
    print '  void Fill(const '+fq_type(klass)+'& sr);'
    print
    print 'protected:'

    # Members
    for tn in basic_mems:
        print '  '+tn.type+' '+tn.name+';'

    if len(basic_mems) > 0 and len(class_mems) > 0: print

    for tn in class_mems:
        print '  '+tn.type+' '+tn.name+';'

    for tn in vec_mems:
        print
        print '  TTree* '+tn.name+'_tree;'
        print '  '+tn.type+' '+tn.name+';'
        print '  long '+tn.name+'_idx;'
        print '  int '+tn.name+'_length;'

    print '};'
    print '} // end namespace'

    # Switch to appending into the main cxx file
    sys.stdout = file(cxxDir+'/FlatRecord.cxx', 'a')

    print
    print '#include "sbncode/FlatMaker/'+pt+'.h"'
    print
    print 'flat::'+pt+'::'+pt+'(const std::string& prefix, TTree* tr, const IBranchPolicy* policy)'

    # Initializer list
    inits = [tn.name+'(prefix+"'+tn.name+'.", tr, policy)' for tn in class_mems]

    for tn in vec_mems:
        inits += [tn.name+'_tree(make_tree(prefix+"'+tn.name+'", "'+tn.name+'", tr))']

        if is_fundamental(tn.type):
            inits += [tn.name+'(0)']
        else:
            inits += [tn.name+'((prefix+"'+tn.name+'."), '+tn.name+'_tree, policy)']

        if tn.is_array():
            inits += [tn.name+'_idx(0), '+tn.name+'_length(sizeof('+fq_type(klass)+'::'+tn.name+')/sizeof(*'+fq_type(klass)+'::'+tn.name+'))']
        else:
            inits += [tn.name+'_idx(0), '+tn.name+'_length(0)']

    if len(inits) > 0: print '  : '+',\n    '.join(inits)

    print '{'

    # Constructor body
    for tn in basic_mems:
        print '  '+branch_cmd(tn.name)

    for tn in vec_mems:
        if is_fundamental(tn.type):
            print '  '+branch_cmd(tn.name, tn.name+'_tree')
        print '  if('+tn.name+'_tree->GetNbranches() > 0){'
        print '    '+branch_cmd_no_policy(tn.name+'_idx')
        print '    '+branch_cmd_no_policy(tn.name+'_length')
        print '  }'

    print '}'

    print

    print 'flat::'+pt+'::~'+pt+'()'
    print '{'
    for tn in vec_mems:
        trName = tn.name+'_tree'
        print '  if('+trName+'->GetNbranches() > 0) '+trName+'->Write();'
        print '  delete '+trName+';'
    print '}'

    print

    print 'void flat::'+pt+'::Fill(const '+fq_type(klass)+'& sr)'
    print '{'
    if klass.name == 'TVector3':
        print '  x = sr.X();'
        print '  y = sr.Y();'
        print '  z = sr.Z();'
        print '}'
        continue

    for tn in basic_mems:
        print '  '+tn.name+' = sr.'+tn.name+';'

    if len(basic_mems) > 0 and len(class_mems) > 0: print

    for tn in class_mems:
        print '  '+tn.name+'.Fill(sr.'+tn.name+');'

    for tn in vec_mems:
        print
        print '  '+tn.name+'_idx += '+tn.name+'_length; // increment taken by previous record'
        if not tn.is_array(): # have to re-fill each time for vectors
            print '  '+tn.name+'_length = sr.'+tn.name+'.size();'
        print '  for(const auto& x: sr.'+tn.name+'){'
        if is_fundamental(tn.type):
            print '    '+tn.name+' = x;'
        else:
            print '    '+tn.name+'.Fill(x);'
        print '    '+tn.name+'_tree->Fill();';
        print '  }'
    print '}'

sys.stderr.write('Wrote FlatRecord.cxx and associated headers\n')
