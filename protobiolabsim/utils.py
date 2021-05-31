"""
Contains various utility methods.
"""



def coalesce ( *arg ) :
    """
    Return the first non-null argument given. Emulates a nullish coalesce operator (PEP 505).
    """
    return next( (a for a in arg if a is not None), None )



def Help_Progressbar(n, loading_time, add):
    '''function for display of a loading bar, n: width of loading bar'''
    import sys
    import time

    loading = '.' * n
    for i in range(n+1):
        # this loop replaces each dot with a hash!
        print('\r%s progress{}: %3d percent'.format(add) % (loading, i*100/n), end='')
        loading = loading[:i] + '#' + loading[i+1:]
        time.sleep(loading_time)
    sys.stdout.write("\n")
