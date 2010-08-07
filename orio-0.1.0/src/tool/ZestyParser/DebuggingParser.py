# ZestyParser 0.8.1 -- Parses in Python zestily
# Copyright (C) 2006-2007 Adam Atlas
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

'''
@version: 0.8.1
@author: Adam Atlas
@copyright: Copyright 2006-2007 Adam Atlas. Released under the MIT license (see LICENSE.txt).
@contact: adam@atlas.st
'''

import Parser, sys

#todo - configurable levels of verbosity?
#use `logging` stdlib module?

class DebuggingParser (Parser.ZestyParser):
    '''
    A L{Parser.ZestyParser} subclass which is useful for debugging parsers. It parses as usual, but it also prints a comprehensive trace to stderr.
    '''
    depth = -1
    
    def __init__(self, *a, **k):
        self.dest = k.pop('dest', sys.stderr)
        Parser.ZestyParser.__init__(self, *a, **k)
    
    def scan(self, token):
        self.depth += 1
        ind = ' |  ' * self.depth
        
        self.dest.write('%sBeginning to scan for %r at position %i\n' % (ind, token, self.cursor))
        r = Parser.ZestyParser.scan(self, token)
        
        if self.last:
            self.dest.write('%sGot %r -- now at %i\n' % (ind, r, self.cursor))
        else:
            self.dest.write("%sDidn't match\n" % ind)
        
        self.depth -= 1
        
        return r
    
    def skip(self, token):
        self.depth += 1
        ind = ' |  ' * self.depth
        
        self.dest.write('%sBeginning to skip %r at position %i\n' % (ind, token, self.cursor))
        r = Parser.ZestyParser.skip(self, token)
        
        if r:
            self.dest.write('%sMatched -- now at %i\n' % (ind, self.cursor))
        else:
            self.dest.write("%sDidn't match\n" % ind)
        
        self.depth -= 1
        
        return r
    
    def iter(self, token, *args, **kwargs):
        self.depth += 1
        ind = ' |  ' * self.depth
        
        self.dest.write('%sBeginning to iterate %r at position %i\n' % (ind, token, self.cursor))
        
        i = Parser.ZestyParser.iter(self, token, *args, **kwargs)
        while 1:
            self.dest.write('%sIterating\n' % ind)
            yield i.next()
        
        self.dest.write('%sDone iterating -- now at %i\n' % (ind, self.cursor))

        self.depth -= 1