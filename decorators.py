
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
    from itertools import starmap

    if kwargs:
        exclude, f = set(kwargs['exclude']), None
        sieve = lambda l:filter(lambda nv: nv[0] not in exclude, l)
    elif len(names) == 1 and isfunction(names[0]):
        f = names[0]
        sieve = lambda l:l
    else:
        names, f = set(names), None
        sieve = lambda l: filter(lambda nv: nv[0] in names, l)
    def decorator(f):
        fargnames, _, _, fdefaults = getargspec(f)
        # Remove self from fargnames and make sure fdefault is a tuple
        fargnames, fdefaults = fargnames[1:], fdefaults or ()
        defaults = list(sieve(zip(reversed(fargnames), reversed(fdefaults))))
        @wraps(f)
        def decorated(self, *args, **kwargs):
            assigned = dict(sieve(zip(fargnames, args)))
            assigned.update(sieve(kwargs.iteritems()))
            for _ in starmap(assigned.setdefault, defaults): pass
            self.__dict__.update(assigned)
            return f(self, *args, **kwargs)
        return decorated
    return f and decorator(f) or decorator


def autoargs(*include,**kwargs):   
    """
    Automatically assigns args and kwargs as attributes
    From http://stackoverflow.com/a/3653049
    """
    import inspect
    import functools
    def _autoargs(func):
        attrs,varargs,varkw,defaults=inspect.getargspec(func)
        def sieve(attr):
            if kwargs and attr in kwargs['exclude']: return False
            if not include or attr in include: return True
            else: return False            
        @functools.wraps(func)
        def wrapper(self,*args,**kwargs):
            # handle default values
            for attr,val in zip(reversed(attrs),reversed(defaults)):
                if sieve(attr): setattr(self, attr, val)
            # handle positional arguments
            positional_attrs=attrs[1:]            
            for attr,val in zip(positional_attrs,args):
                if sieve(attr): setattr(self, attr, val)
            # handle varargs
            if varargs:
                remaining_args=args[len(positional_attrs):]
                if sieve(varargs): setattr(self, varargs, remaining_args)                
            # handle varkw
            if kwargs:
                for attr,val in kwargs.iteritems():
                    if sieve(attr): setattr(self,attr,val)            
            return func(self,*args,**kwargs)
        return wrapper
    return _autoargs
