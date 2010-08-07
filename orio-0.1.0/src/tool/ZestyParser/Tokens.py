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
@group Basic Tokens: Raw,RE,RawToken,Token,TakeToken
@group Complex Tokens: CompositeToken,TokenSequence,TokenSeries
@group Special Tokens: Defer,Default,Lookahead,Negative
@group TokenSequence Flags: Omit,Skip,Only

@version: 0.8.1
@author: Adam Atlas
@copyright: Copyright 2006-2007 Adam Atlas. Released under the MIT license (see LICENSE.txt).
@contact: adam@atlas.st

@var EmptyToken: A L{Default} instance initialized with the empty string.
@var EOF: A token which matches (and returns C{None}) if the parser is at the end of its L{data <ZestyParser.data>} sequence.

In ZestyParser, a token object must, at minimum, be a callable taking a L{ZestyParser <Parser.ZestyParser>} instance and its current L{cursor <ZestyParser.cursor>} as parameters. It can do whatever it needs with the parser's L{data <ZestyParser.data>} and L{cursor <ZestyParser.cursor>} properties before returning. It may raise L{NotMatched} to indicate to the L{ZestyParser <Parser.ZestyParser>} instance that it failed to match; it may also raise L{ParseError} to indicate, for instance, that it began matching successfully but encountered an unrecoverable error.

The L{Tokens} module contains a variety of predefined token classes (instances of which are callable) and other valid token objects which should cover most parsing situations.
'''

import re, copy, types, warnings
from Parser import NotMatched, ParseError

__all__ = ('Placeholder', 'AbstractToken', 'TokenWrapper', 'RE', 'Raw', 'Token', 'RawToken', 'CompositeToken', 'TokenSequence', 'TakeToken', 'TokenSeries', 'EmptyToken', 'Default', 'Skip', 'Omit', 'Only', 'Defer', 'Lookahead', 'Negative', 'EOF', 'Whitespace', 'Const', 'Inf')


rstack = []
replstack = []

Inf = -1

def count_args(callable):
    t = type(callable)
    if t is types.FunctionType:
        return callable.func_code.co_argcount
    elif t is types.ClassType:
        return callable.__init__.im_func.func_code.co_argcount - 1
    elif t is types.InstanceType:
        return callable.__call__.im_func.func_code.co_argcount - 1
    elif t is types.MethodType:
        return callable.im_func.func_code.co_argcount - 1
    #assume it's some builtin that only takes the data itself as a parameter
    return 1

class Placeholder:
    def __init__(self, key=None):
        self.key = key
    
    def __eq__(self, key):
        return self.key == key
    
    def __hash__(self):
        return hash(self.key)
    
    def __call__(self, key=None):
        return Placeholder(key)
    
    def _single_replace(cls, subj, vals, kwvals):
        if isinstance(subj, cls):
            if vals and subj == None:
                return vals.pop(0)
            elif subj in kwvals and subj != None:
                return kwvals[subj]
        elif isinstance(subj, AbstractToken):
            subj._replace(vals, kwvals)
        return subj
    _single_replace = classmethod(_single_replace)
    
    def _list_replace(cls, subj, vals, kwvals):
        for i, v in enumerate(subj):
            subj[i] = cls._single_replace(v, vals, kwvals)
    _list_replace = classmethod(_list_replace)
    
    def __repr__(self):
        return '<Placeholder %r>' % self.key

class ListReplacing:    
    def _replace(self, vals, kwvals):
        if self not in replstack:
            replstack.append(self)
            Placeholder._list_replace(self.desc, vals, kwvals)
            replstack.pop()
class SingleReplacing:        
    def _replace(self, vals, kwvals):
        if self not in replstack:
            replstack.append(self)
            self.desc = Placeholder._single_replace(self.desc, vals, kwvals)
            replstack.pop()

class AbstractToken (object):
    '''
    Base class from which most tokens defined in this module derive. Subclassing this is not required for writing tokens, since they can be any callable with certain semantics, but this class provides several useful services for creating reusable token classes, such as callback support and convenient operator overloading.
    
    @ivar desc: The generic "description" variable which stores the "essence" of any given instance. Subclasses use this as needed.
    @ivar callback: An optional callable which, if not None, will be called whenever an instance matches successfully. It may take one, two, or three parameters, depending on its needs. If one, it will be passed whatever data the token matched (i.e. whatever it would normally have returned upon being called). If two, it will be passed the L{ZestyParser <Parser.ZestyParser>} instance and the data. If three, it will be passed the parser, the data, and the what the parser's cursor was when this token started matching. Callbacks may raise L{NotMatched} or L{ParseError} with the usual behaviour. They should also return a value, which will be returned to the calling L{ZestyParser <Parser.ZestyParser>} instance.
    @ivar to: An optional callable which, if not None, will be called in the same manner as a callback (after any callback and before returning to the parser instance), but will be passed only one argument: the data matched (or returned by the callback, if any). Its main purpose is to allow you to concisely do things like C{Token('[0-9]+', group=0, to=int)} -- the builtin callable C{int} will be passed the text matched by the regex, so the token will ultimately return an integer instead of a string or a regex match object. You can also use this property with L{AHT} types, for more complex multi-stage parsing. See the C{n3.py} and C{n3rdflib.py} examples for a demonstration of this. (In previous versions, this was passed to the initializer as C{as}, but this is deprecated because C{as} will become a reserved word in Python 2.6. Change your code to use {to}.)
    '''
    
    name = None
    failMessage = None
    callback = None
    to = None

    #'as' is deprecated in favor of 'to' since it's becoming a reserved word
    def __init__(self, desc, callback=None, to=None, name=None):
        self.desc = desc
        self.callback = callback
        self.to = to 
        self.name = name
    
    def __repr__(self):
        return '%s %s' % (self.__class__.__name__, (self.name or str(self)))

    def __str__(self): return ''

    def _make_callbackrun(self, func, callback):
        argcount = count_args(callback)
        if argcount == 1:
            def f(parser, origCursor):
                return callback(func(parser, origCursor))
        elif argcount == 2:
            def f(parser, origCursor):
                return callback(parser, func(parser, origCursor))
        elif argcount == 3:
            def f(parser, origCursor):
                return callback(parser, func(parser, origCursor), origCursor)
        return f
    
    def _make_torun(self, func):
        def f(parser, origCursor):
            return self.to(func(parser, origCursor))
        return f
    
    def _make_failcheck(self, func):
        def f(parser, origCursor):
            try:
                data = func(parser, origCursor)
                return data
            except NotMatched:
                if self.failMessage is True:
                    raise ParseError(parser, 'Expected %s' % str(self))
                elif self.failMessage:
                    raise ParseError(parser, self.failMessage)
                else: raise
        return f
    
    def _poke(self):
        c = self.__call__
        if self.callback:
            c = self._make_callbackrun(c, self.callback)
        if self.to:
            c = self._make_torun(c)
        if self.failMessage:
            c = self._make_failcheck(c)
        
        if c is self.__call__ and not isinstance(self, Defer):
            self.parse = None
            del self.parse
        else:
            self.parse = c
        
    def __copy__(self):
        n = self.__class__.__new__(self.__class__)
        n.__dict__.update(self.__dict__)
        n._poke()
        n.desc = copy.copy(self.desc)
        return n
    
    def __setattr__(self, name, value):
        super(AbstractToken, self).__setattr__(name, value)
        if name in ('callback', 'failMessage', 'to'):
            self._poke()

    def __add__(self, other):
        '''Allows you to construct L{TokenSequence}s with the + operator.'''
        return TokenSequence([self, other])
    
    def __sub__(self, other):
        '''Allows you to construct L{TokenSequence}s with the - operator, automatically padded with L{Whitespace}.
        
I realize it's a bit weird to use the - operator for this, but the main motivation is giving it the same precedence as +. Still, you can read it as a sort of "blank" (which is what the left and right tokens are being joined by), instead of "minus".'''
        return TokenSequence([self, Whitespace, other])
    
    def __or__(self, other):
        '''Allows you to construct L{CompositeToken}s with the | operator.'''
        return CompositeToken([self, other])
    
    def __mul__(self, val):
        '''Allows you to construct L{TokenSeries} with the * operator. Operand can be:
        
            - int (a series of exactly this many)
            - (int, ) (a series of at least this many)
            - (x:int, y:int) a series of x to y
        
        The constant Inf can be used in some of these -- * Inf yields a 0--infinity series, and * (x, Inf) yields an x--infinity series.
        '''
        t = TokenSeries(self)
        if isinstance(val, int):
            t.min = t.max = val
        elif isinstance(val, tuple):
            if len(val) == 2:
                t.min, t.max = val
            elif len(val) == 1:
                t.min = val
        return t
    __rmul__ = __mul__
    
    def __rshift__(self, callback):
        '''
        Convenience overloading for setting the L{callback<AbstractToken.callback>} of a token whose initializer you do not call directly, such as the result of combining tokens with L{+<__add__>} or L{|<__or__>}.
        
        @param callback: An L{AbstractToken}-compatible callback.
        @type callback: callable
        @return: A copy of C{self} with the L{callback<AbstractToken.callback>} property set to C{callback}.
        '''
        new = copy.copy(self)
        new.callback = callback
        return new
    
    def __xor__(self, message):
        '''
        Overloading for setting the L{failMessage<AbstractToken.failMessage>} of a token.
        
        @param message: The message to be raised with L{ParseError} if this token fails to match.
        @type message: str
        @return: A copy of C{self} with the L{failMessage<AbstractToken.failMessage>} property set to C{callback}.
        '''
        new = copy.copy(self)
        new.failMessage = message
        return new
    
    def __invert__(self):
        return Negative(self)
        
    def _replace(self, vals, kwvals):
        pass
    
    def __imod__(self, val):
        if isinstance(val, (tuple, list)):
            self._replace(val, {})
        elif isinstance(val, dict):
            self._replace([], val)
        else:
            self._replace([val], {})
        return self
    
    def __mod__(self, val):
        new = copy.copy(self)
        new %= val
        return new

class TokenWrapper (AbstractToken):
    '''If you write your own token type in a way other than subclassing AbstractToken, e.g. by simply writing a function, you can use this as a decorator to automatically let it take advantage of AbstractToken's magic.
    '''
    def __call__(self, parser, origCursor):
        t = parser.scan(self.desc)
        if parser.last:
            return t
        else: raise NotMatched

    def __str__(self):
        return repr(self.desc)

class RE (AbstractToken):
    '''
    A class whose instances match Python regular expressions.
    
    @ivar group: If defined, L{__call__} returns that group of the regular expression match instead of the whole match object.
    @type group: int
    '''
    
    def __init__(self, regex, group=None, **kwargs):
        '''
        @param regex: Either a compiled regex object or a string regex.
        @param group: To be set as the object's L{group} property.
        @type group: int
        '''
        if not hasattr(regex, 'match'):
            regex = re.compile(regex, re.DOTALL)
        super(Token, self).__init__(regex, **kwargs)
        if group is not None:
            try:
                group = int(group)
            except ValueError:
                raise ValueError('got non-numeric value for `group`: ' + group)
        self.group = group
    
    def __call__(self, parser, origCursor):
        matches = self.desc.match(parser.data, origCursor)
        if matches is None: raise NotMatched

        parser.cursor = matches.end()
        if self.group is not None:
            matches = matches.group(self.group)
        return matches
    
    def __str__(self):
        return repr(self.desc.pattern)

Token = RE

class Raw (AbstractToken):
    '''
    A class whose instances match only a particular string. Returns that string.
    
    @ivar caseInsensitive: If true, ignores case.
    @type caseInsensitive: bool
    '''
    def __init__(self, string, caseInsensitive=False, **kwargs):
        '''
        @param string: The string to match.
        @type string: str
        @param caseInsensitive: To be set as the object's L{caseInsensitive} property.
        @type caseInsensitive: bool
        '''
        super(Raw, self).__init__(string, **kwargs)
        self.len = len(string)
        self.caseInsensitive = caseInsensitive
        if caseInsensitive:
            self.desc = self.desc.lower()
    
    def __call__(self, parser, origCursor):
        end = origCursor + self.len
        d = parser.data[origCursor:end]
        if (not self.caseInsensitive and d == self.desc) or (self.caseInsensitive and d.lower() == self.desc):
            parser.cursor = end
            return d
        else: raise NotMatched
    
    def __str__(self):
        return repr(self.desc)

RawToken = Raw

class Default (AbstractToken):
    '''
    A class whose instances always return L{desc} and do not advance the parser's cursor.
    '''
    def __call__(self, parser, origCursor):
        return self.desc
    def __str__(self):
        return repr(self.desc)

EmptyToken = Default('')

class CompositeToken (ListReplacing, AbstractToken):
    '''
    A class whose instances match any of a number of tokens.
    
    @ivar desc: A list of token objects.
    @type desc: list
    '''
    def __call__(self, parser, origCursor):
        for t in self.desc:
            r = parser.scan(t)
            if parser.last:
                return r
        raise NotMatched
    
    def __str__(self):
        if self in rstack:
            return '...'
        else:
            rstack.append(self)
            d = '(' + ' | '.join([repr(t) for t in self.desc]) + ')'
            rstack.pop()
            return d
    
    def __or__(self, other):
        if hasattr(other, '__iter__'):
            return CompositeToken(self.desc + list(other))
        else:
            return CompositeToken(self.desc + [other])
    
    def __ior__(self, other):
        if hasattr(other, '__iter__'):
            self.desc += list(other)
        else:
            self.desc.append(other)
        return self

class TokenSequence (ListReplacing, AbstractToken):
    '''
    A class whose instances match a sequence of tokens. Returns a corresponding list of return values from L{ZestyParser.scan}.
    
    Some special types, L{Skip}, L{Omit}, and L{Only}, are allowed in the sequence. These are wrappers for other token objects adding special behaviours. If it encounters a L{Skip} token, it will process it with L{ZestyParser.skip}, ignore whether it matched, and not include it in the list. If it encounters a L{Omit} token, it will still require that it match (the default behaviour), but it will not be included in the list.
    
    If the sequence contains an L{Only} token, its result will be returned instead of the usual list, though it still requires that subsequent tokens match. Multiple L{Only} tokens are meaningless and L{TokenSequence}'s behavior in that case is undefined. 
    
    @ivar desc: A list of token objects.
    @type desc: list
    '''
    
    def __call__(self, parser, origCursor):
        o = []
        only = False
        onlyVal = None
        for g in self.desc:
            if g is Whitespace:
                parser.skip(parser.whitespace)
            r = parser.scan(g)
            if parser.last is None:
                raise NotMatched
            if isinstance(g, Only):
                only = True
                onlyVal = r
                continue
            if not isinstance(g, (Skip, Omit, _Whitespace)) and not only:
                o.append(r)
        if only: #heh
            return onlyVal
        else:
            return o
        
    def __str__(self):
        if self in rstack:
            return '...'
        else:
            rstack.append(self)
            d = '(' + ' + '.join([repr(t) for t in self.desc]) + ')'
            rstack.pop()
            return d

    def __add__(self, other):
        if hasattr(other, '__iter__'):
            return TokenSequence(self.desc + list(other))
        else:
            return TokenSequence(self.desc + [other])
    
    def __sub__(self, other):
        if hasattr(other, '__iter__'):
            return TokenSequence(self.desc + [Whitespace] + list(other))
        else:
            return TokenSequence(self.desc + [Whitespace, other])
    
    def __iadd__(self, other):
        if hasattr(other, '__iter__'):
            self.desc += list(other)
        else:
            self.desc.append(other)
        return self

class TakeToken (AbstractToken):
    '''
    A class whose instances match and return a given number of characters from the parser's L{data<ZestyParser.data>}. Raises L{NotMatched} if not enough characters are left.
    '''
    
    def __init__(self, length, **kwargs):
        super(TakeToken, self).__init__(length, **kwargs)
    
    def __call__(self, parser, start):
        end = start + self.desc
        if parser.len < end: raise NotMatched
        parser.cursor = end
        return parser.data[start:end]

class TokenSeries (SingleReplacing, AbstractToken):
    '''
    A particularly versatile class whose instances match one token multiple times (with a great degree of customizability).
    
    The properties L{skip}, L{prefix}, L{postfix}, and L{delimiter} are optional tokens which add structure to the series. It can be represented, approximately in the idioms of L{TokenSequence}, as follows::
    
        [Skip(skip) + Omit(prefix) + desc + Omit(postfix)] + [Skip(skip) + Omit(delimiter) + Skip(skip) + Omit(prefix) + desc + Omit(postfix)] + ... + Skip(skip)
    
    Or, if there is no delimiter::
    
        [Skip(skip) + Omit(prefix) + desc + Omit(postfix)] + ... + Skip(skip)
    
    @ivar desc: The token to match.
    @type desc: token
    @ivar min: The minimum number of times L{desc} must match.
    @type min: int
    @ivar max: The maximum number of times L{desc} will try to match.
    @type max: int
    @ivar skip: An optional token to skip between matches.
    @type skip: token
    @ivar prefix: An optional token to require (but omit from the return value) before each instance of L{token}.
    @type prefix: token
    @ivar postfix: An optional token to require (but omit from the return value) after each instance of L{token}.
    @type postfix: token
    @ivar delimiter: An optional token to require (but omit from the return value) between each instance of L{token}.
    @type delimiter: token
    @ivar until: An optional 2-tuple whose first item is a token, and whose second item is either a message or False. The presence of this property indicates that the token in C{until[0]} must match at the end of the series. If this fails, then if C{until[1]} is a message, a ParseError will be raised with that message; if it is False, NotMatched will be raised.
    '''
    def __init__(self, token, min=0, max=-1, skip=None, prefix=None, postfix=None, delimiter=None, until=None, includeDelimiter=False, **kwargs):
        super(TokenSeries, self).__init__(token, **kwargs)
        self.min, self.max, self.skip, self.prefix, self.postfix, self.delimiter, self.until, self.includeDelimiter = min, max, skip, prefix, postfix, delimiter, until, includeDelimiter

    def __call__(self, parser, origCursor):
        o = []
        i = 0
        done = False
        while i != self.max:
            c = parser.cursor
            if self.until and parser.skip(self.until[0]):
                done = True
                break
            if self.skip: parser.skip(self.skip)

            c = parser.cursor
            if i != 0 and self.delimiter is not None:
                d = parser.scan(self.delimiter)
                if parser.last is None: break
                if self.skip: parser.skip(self.skip)
            if self.prefix and not parser.skip(self.prefix): break
            t = parser.scan(self.desc)
            if parser.last is None: break
            if self.postfix and not parser.skip(self.postfix): break
            
            if i != 0 and self.includeDelimiter: o.append(d)
            o.append(t)
            i += 1
        if self.until and not done:
                parser.cursor = c
                if self.until[1]:
                    if self.until[1] is True:
                        raise ParseError(parser, 'Expected %s' % self.until[0])
                    else:
                        raise ParseError(parser, self.until[1])
                else:
                    raise NotMatched
        if i >= self.min:
            return o
        else:
            raise NotMatched
    
    def __str__(self):
        return repr(self.desc)

class Defer (AbstractToken):
    '''
    A token which takes a callable (generally a lambda) which takes no arguments and itself returns a token. A Defer instance, upon being called, will call this function, scan for the returned token, and return that return value. This is primarily intended to allow you to define tokens recursively; if you need to refer to a token that hasn't been defined yet, simply use C{Defer(lambda: T_SOME_TOKEN)}, where C{T_SOME_TOKEN} is the token's eventual name.
    '''
    
    def __call__(self, parser, origCursor):
        self.__call__ = self.desc()
        self._poke()
        return self.parse(parser, origCursor)

def _pad(p, subj, outer=True):
    if isinstance(subj, TokenSequence):
        t = copy.copy(subj)
        t.desc = list(t.desc)
        start = (not outer) and 1 or 0
        stop = 2*(len(t.desc) + (outer and 1 or -1))
        for i in range(start, stop, 2): t.desc.insert(i, p)
        return t
    else:
        return TokenSequence([p, Only(subj), p])

class Skip (SingleReplacing, AbstractToken):
    '''
    See L{TokenSequence}.
    '''
    def __call__(self, parser, origCursor):
        parser.skip(self.desc)
    
    def pad(self, tok, outer=True):
        '''
        Takes a TokenSequence and returns a copy with this L{Skip} token separating every element therein. Alternately, takes any other token and returns a TokenSequence equivalent to (self + Only(tok) + self).
        
        @param tok: The token sequence to operate on.
        @type tok: L{TokenSequence} or other token
        @param outer: If operating on a L{TokenSequence}, whether to also include self on either side of the sequence, in addition to between each element.
        @type outer: bool
        @return: L{TokenSequence}
        '''
        return _pad(self, tok, outer)

class Omit (SingleReplacing, AbstractToken):
    '''
    See L{TokenSequence}.
    '''
    def __call__(self, parser, origCursor):
        if not parser.skip(self.desc):
            raise NotMatched

class Only (SingleReplacing, TokenWrapper):
    '''
    See L{TokenSequence}.
    '''

class Lookahead (SingleReplacing, AbstractToken):
    '''
    Scans for another token and returns its result as usual, but doesn't actually advance the parser's cursor.
    '''
    def __call__(self, parser, origCursor):
        t = parser.scan(self.desc)
        parser.cursor = origCursor
        if parser.last is None: raise NotMatched
        return t

class Negative (SingleReplacing, AbstractToken):
    '''
    Scans for another token, and only matches (returning True) if that token did not match.
    '''
    def __call__(self, parser, origCursor):
        parser.scan(self.desc)
        parser.cursor = origCursor
        if parser.last: raise NotMatched
        return True

class _EOF (AbstractToken):
    '''
    Matches returning None if the parser is at the end of its input.
    '''
    def __call__(self, parser, origCursor):
        if parser.cursor != parser.len: raise NotMatched
EOF = _EOF(None)

class _Whitespace (AbstractToken):
    def __call__(self, parser, origCursor):
        parser.skip(parser.whitespace)
Whitespace = _Whitespace(None)

def Const(value):
    '''
    This is a factory function used for when you want an L{AbstractToken}-compatible callback that returns a constant. Example::
        RawToken('foo') >> Const('bar')
    
    This token matches 'foo' as usual, but always returns 'bar'.
    '''
    def f(r):
        return value
    return f
