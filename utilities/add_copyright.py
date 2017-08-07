#!/usr/bin/env python
import fileinput;
import sys;
import re;
import os;


def enum(*members):
    """Emulation of enums

    Usage:
    e=enum(ZERO,ONE,TWO,THREE)
    assert e.ONE==1
    """
    return type('ENUM', (), dict(zip(members, range(0,len(members)))))


def process(source,print_licence,is_empty,is_begin_mark,is_end_mark):
    """Read from FileInput source, insert licence between markers
    
    Arguments:
    source --  FileInput source (must be iterator and support ``isfirstline()`` method)
    print_licence -- a callable that prints the licence
    is_empty -- callable predicate detecting "empty" lines that may precede the licence
    is_begin_mark -- callable predicate detecting "licence begins" marker
    is_begin_mark -- callable predicate detecting "licence ends" marker
    """
    # List the states of the parser
    States=enum('SEARCHING_HEADER', 'SKIPPING_TO_FOOTER','AFTER_FOOTER')

    state=States.SEARCHING_HEADER
    try:
        line=source.next()
        while True: # loop will be exited on exception
            if state==States.SEARCHING_HEADER:
                if is_begin_mark(line):
                    state=States.SKIPPING_TO_FOOTER
                    line=source.next()
                    continue
                if is_empty(line):
                    print line,
                    line=source.next()
                    continue
                print_licence()
                print line,
                state=States.AFTER_FOOTER
                continue
            if state==States.SKIPPING_TO_FOOTER:
                if source.isfirstline():
                    state=States.SEARCHING_HEADER
                    continue
                if is_end_mark(line):
                    print_licence()
                    state=States.AFTER_FOOTER
                    continue
                line=source.next()
                continue
            if state==States.AFTER_FOOTER:
                line=source.next()
                if source.isfirstline():
                    state=States.SEARCHING_HEADER
                else:
                    print line,
                continue
            assert False,"Impossible state #"+str(state)

    except StopIteration:
        pass
    return


def Main():
    if len(sys.argv)<3:
        sys.stderr.write("""Add copyright statement to all *.cpp and *.hpp files in the given directory trees.

Usage: %s copyright_file dirtree [dirtree...]
""" % sys.argv[0])
        return 2

    (cpr_name, dir_names)=(sys.argv[1], sys.argv[2:])

    # Sentries for the copyright statement:
    cpr_start="/*** LICENCE: ***"
    cpr_end="*** END OF LICENCE ***/"
    
    # Prepare copyright blurb
    try:
        cpr_file=open(cpr_name,"r")
    except IOError:
        sys.stderr.write("Cannot open file '%s'\n" % cpr_name)
        return 1
    cpr_text=cpr_start+"\n"+cpr_file.read()+"\n"+cpr_end+"\n"
    cpr_file.close()

    # make inner functions to print the copyright text, catch start/stop, etc
    def print_licence():
        print cpr_text,

    def is_empty(line):
        return re.search(r"^\s*$", line)!=None
    
    def is_begin_mark(line):
        return re.search(r"^\s*"+re.escape(cpr_start), line)!=None
    
    def is_end_mark(line):
        return re.search(re.escape(cpr_end), line)!=None

    # Get the list of files to be processed
    files=[]
    for topdir in dir_names:
        for base_dirs_fnames in os.walk(topdir):
            for fname in base_dirs_fnames[2]:
                if fname[-4:]=='.cpp' or fname[-4:]=='.hpp':
                    files.append(os.path.join(base_dirs_fnames[0], fname))
    
    # Process the files, replacing their content
    process(fileinput.FileInput(files,inplace=True), print_licence, is_empty, is_begin_mark, is_end_mark)

    return 0

# Run the Main and exit
sys.exit(Main())
