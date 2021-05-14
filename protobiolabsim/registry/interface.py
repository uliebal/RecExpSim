"""
Registries provide an efficient way to store Records that are used by multiple Organisms.
That way, each organism only needs to store a reference to each registered object (for example
a gene) instead to storing and copying the entire gene.
When editing an object in the registry, the editor should make sure to practice copy-on-write.

Entries of a registry must be copy-able.
"""

from abc import ABC, abstractmethod
from typing import TypeVar, Generic, List
from copy import copy



T = TypeVar('T')



# class T ( ABC ) : #Generic[T],
#     """ T is something that is stored in a Registry. """
#
#     @abstractmethod
#     def __copy__ ( self ) :
#         pass



class Registry ( Generic[T] ) :
    """ A Registry holds multiple Entries. """

    records : set[T]

    def __init__ ( self ) :
        self.records = set()


    def insert ( self, entry: T ) -> None :
        self.records.add(entry)


    def clone ( self, entry: T ) -> T :
        if entry in self.records :
            copied_entry = copy(entry)
            self.insert(copied_entry)
            return copied_entry
        raise Exception("Tried to clone a non-registered entry.")

