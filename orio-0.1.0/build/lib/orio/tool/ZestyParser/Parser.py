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

__all__ = ('ZestyParser', 'NotMatched', 'ParseError', 'CallbackFor', 'DebugNote')

class Error (Exception): pass
class NotMatched (Error): '''Raised by a token if it has failed to match at the parser's current cursor.'''
class ParseError (Error):
    '''Raised by a token to indicate that a parse error has occurred.'''
    def __init__(self, parser, message):
        '''
        @param parser: The parser instance that encountered the error.
        @type parser: ZestyParser
        @param message: A message explaining the error.
        @type message: str
        '''
        self.parser, self.message, self.coord = parser, message, parser.coord()
    def __str__(self):
        '''Prints the error message and the row and column corresponding to the parser's cursor.'''
        return "%s at line %i column %i" % (self.message, self.coord[0], self.coord[1])

def CallbackFor(token):
    '''
    Function decorator indicating that the function should be set as the callback of the given token; returns the token instead of the function.
    
    Example::
        @CallbackFor(Token('([0-9]+)'))
        def T_INT(r):
            print r
    
    This is equivalent to::
        def T_INT(r):
            print r
        T_INT = Token('([0-9]+)', callback=T_INT)
    '''
    def newFunc(func):
        return token >> func
    return newFunc

def DebugNote(note):
    '''
    Factory which creates a logging function usable as a token callback. The logging function prints L{note}, the parser's coordinates, and the value matched. It then returns the value unmodified.
    
    Example::
        T_FOO = TokenSeries(RawToken('foo')) >> DebugNote('foo_series')
    '''
    def func(parser, val, cursor):
        coord = parser.coord()
        print '%s(%i, %i): %s' % (note, coord[0], coord[1], val)
        return val
    return func

class ZestyParser:
    '''
    Parses one stream of data, by means of L{tokens<ZestyParser.Tokens>}.
    
    @ivar context: A dictionary which can be used for storing any necessary state information.
    @type context: dict
    @ivar data: The sequence being parsed (probably a string).
    @type data: sequence
    @ivar cursor: The current position of the parser in L{data}.
    @type cursor: int
    @ivar last: The last matched token.
    @type last: L{token<ZestyParser.Tokens>}
    '''
    
    context = {}
    data = None
    cursor = 0
    len = 0
    last = None
    
    whitespace = None

    def __init__(self, data=None):
        '''Initializes the parser, optionally calling L{useData}'''
        if data: self.useData(data)
        from Tokens import RE
        self.whitespace = RE('\s+')

    def useData(self, data):
        '''
        Begin parsing a stream of data
        
        @param data: The data to parse.
        @type data: sequence
        '''
        self.data = data
        self.cursor = 0
        self.len = len(data)
    
    def scan(self, token):
        '''
        Scan for one token.
        
        @param token: The token to scan for.
        @return: The return value of the matching token, or None if the token raised NotMatched.
        @rtype: object
        @raise ParseError: If a token fails to match and it has a failMessage parameter.
        '''
        oldCursor = self.cursor
        try:
            r = getattr(token, 'parse', token)(self, oldCursor)
            self.last = token
            return r
        except NotMatched:
            self.cursor = oldCursor
            self.last = None
            return None
    
    def skip(self, token):
        '''
        A convenience method that skips one token and returns whether it matched.
        
        @param token: The token to scan for.
        @type token: token
        @return: Whether or not the token matched.
        @rtype: bool
        '''
        oldCursor = self.cursor
        try:
            getattr(token, 'parse', token)(self, oldCursor)
            return token
        except NotMatched:
            self.cursor = oldCursor
    
    def iter(self, token, skip=None, until=None):
        '''
        Returns a generator iterator which scans for L{token} every time it is invoked.
        
        @param token: The token to scan for.
        @param skip: An optional token to L{skip} before each L{scan} for L{token}.
        @type skip: token
        @param until: An optional 2-tuple. If defined, the iterator will scan for L{token} until it reaches the token C{until[0]}; if L{scan} returns C{None} before the iterator encounters this token, it raises a L{ParseError} with the message given in C{until[1]}.
        @type until: tuple
        @rtype: iterator
        '''
        while 1:
            if skip: self.skip(skip)
            if until and self.skip(until[0]): break
            r = self.scan(token)
            if self.last: yield r
            elif until:
                if until[1] is True:
                    raise ParseError(self, 'Expected %s' % str(until[0]))
                else:
                    raise ParseError(self, until[1])
            else: break
    
    def coord(self, loc=None):
        '''
        Returns row/column coordinates for a given point in the input stream, or L{cursor} by default. Counting starts at C{(1, 1)}.
        
        @param loc: An index of L{data}.
        @type loc: int
        @return: A 2-tuple representing (row, column).
        @rtype: tuple
        '''
        if loc is None: loc = self.cursor
        row = self.data.count('\n', 0, loc) + 1
        col = loc - self.data.rfind('\n', 0, loc)
        return (row, col)