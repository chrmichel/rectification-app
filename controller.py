'''
Main program for the app
version 0.4

23.02.2023

Christian Michel
'''

from typing import Protocol
from configparser import ConfigParser

from rectification import Rectification


class UI(Protocol):
    '''
    Protocol class to define UI specifications.
    '''

    def get_params(self):
        ...
    def use_ini(self):
        ...
    def calc(self, params):
        ...
    def reset(self):
        ...
    def output(self):
        ...


class Controller:
    '''
    Controller connects UI and calculations.
    '''

    def __init__(self, ui: UI) -> None:
        self.ui = ui
        self.params = None
        self.rect = None

    def _read_config(self) -> None:
        cfg = ConfigParser()
        cfg.read("config.ini")
        input_dict = dict(cfg.items("DATA"))
        input_dict["comps"] = [i.strip() for i in input_dict["comps"].split()]
        input_dict["z_f"] = [float(i.strip())
                                for i in input_dict["z_f"].split()]
        self.params= input_dict
        self.ui.use_ini()
        if self.params:
            self.rect = Rectification(self.params)
        else:
            self._reset()

    def _reset(self) -> None:
        self.ui.reset()

    def _get_params(self):
        self.params = self.ui.get_params()
        if self.params:
            self.rect = Rectification(self.params)
        else:
            self._reset()

    def _calc(self) -> None:
        self.ui.calc(self.params)
        if self.rect:
            self.rect.run()

    def _output(self) -> None:
        if self.rect:
            self.ui.output()
