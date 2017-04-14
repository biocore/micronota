from abc import ABC, abstractmethod
from os.path import basename, splitext, join

class BaseMod(ABC):
    def __init__(self, directory, name):
        self.name = join(directory, splitext(basename(name))[0])
        self.result = {}
        self.report = {}

    @abstractmethod
    def parse(self):
        '''parse result'''
