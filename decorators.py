
def autoassign(*names, **kwargs):
    """Decorator that automatically assigns keywords as atributes
 
    allow a method to assign (some of) its arguments as attributes of
    'self' automatically.  E.g.
    
    To restrict autoassignment to 'bar' and 'baz', write:
    
    @autoassign('bar', 'baz')
    def method(self, foo, bar, baz): ...

    To prevent 'foo' and 'baz' from being autoassigned, use:

    @autoassign(exclude=('foo', 'baz'))
    def method(self, foo, bar, baz): ...
    """
    from functools import wraps
    from inspect import getargspec, isfunction
    from itertools import izip, ifilter, starmap

    if kwargs:
        exclude, f = set(kwargs['exclude']), None
        sieve = lambda l:ifilter(lambda nv: nv[0] not in exclude, l)
    elif len(names) == 1 and isfunction(names[0]):
        f = names[0]
        sieve = lambda l:l
    else:
        names, f = set(names), None
        sieve = lambda l: ifilter(lambda nv: nv[0] in names, l)
    def decorator(f):
        fargnames, _, _, fdefaults = getargspec(f)
        # Remove self from fargnames and make sure fdefault is a tuple
        fargnames, fdefaults = fargnames[1:], fdefaults or ()
        defaults = list(sieve(izip(reversed(fargnames), reversed(fdefaults))))
        @wraps(f)
        def decorated(self, *args, **kwargs):
            assigned = dict(sieve(izip(fargnames, args)))
            assigned.update(sieve(kwargs.iteritems()))
            for _ in starmap(assigned.setdefault, defaults): pass
            self.__dict__.update(assigned)
            return f(self, *args, **kwargs)
        return decorated
    return f and decorator(f) or decorator

