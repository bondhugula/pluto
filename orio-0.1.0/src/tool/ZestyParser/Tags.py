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

Tags is a utility module providing an easy way to label objects in abstract parse trees without defining a class for each one. It supersedes the L{AHT} module.

This module provides a global "Tags" object; you create a tag by accessing any attribute on it. A tag is a callable object, suitable for use as an L{AbstractToken} C{to} parameter, or, if it is more convenient (e.g. when you must use C{>>}), a callback. Later, you can check if a given tag has been applied to an object by checking for membership with C{in}. For example::

    >>> l = [1, 2, 3]
    >>> l in Tags.thing
    False
    >>> Tags.thing(l)
    [1, 2, 3]
    >>> l in Tags.thing
    True
'''

__all__ = ('Tags',)

class _Env (object):
    '''
    @see: L{Tags}
    '''
    _tagobjs = {}
    
    def __getattr__(self, attr):
        if attr not in self._tagobjs:
            self._tagobjs[attr] = _Tag(attr)
        return self._tagobjs[attr]
Tags = _Env()

class _Tag (object):
    def __init__(self, name):
        self.name = name
        self.objects = []
    
    def __call__(self, arg):
        self.objects.append(arg)
        return arg
    
    def __contains__(self, arg):
        return (arg in self.objects)
    
    def __repr__(self):
        return '<Tag %s>' % self.name