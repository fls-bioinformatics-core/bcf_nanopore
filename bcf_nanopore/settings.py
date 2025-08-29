#!/bin/env python
#
#     settings.py: handle package configuration
#     Copyright (C) University of Manchester 2025 Peter Briggs
#

"""
Module for handling package configuration information from
settings for automated processing from '.ini'-style config
files.

Provides the following classes:

- Settings: handle configuration for 'bcf_nanopore' package
"""

#######################################################################
# Imports
#######################################################################

import os
from auto_process_ngs.settings import GenericSettings
from auto_process_ngs.settings import locate_settings_file
from auto_process_ngs.settings import get_config_dir
from auto_process_ngs.settings import jobrunner

#######################################################################
# Classes
#######################################################################

class Settings(GenericSettings):
    """
    Handle local configuration for ``bcf_nanopore``

    By default tries to locate a configuration file
    called ``bcf_nanopore.ini`` by checking the
    following locations:

    1. Current directory
    2. Config directory of the installation

    Arguments:
      settings_file (str): path to configuration file to
        read settings from (optional)
      resolve_undefined (bool): if True (default) then
        assign values to "null" parameters by checking
        fallback parameters and default values
    """
    def __init__(self, settings_file=None, resolve_undefined=True):
        GenericSettings.__init__(
            self,
            # Define the sections, parameters and types
            settings = {
                "general": { "default_runner": jobrunner,
                             "permissions": str,
                             "group": str },
                "runners": { "rsync": jobrunner },
                "reporting_templates": { "*": str },
            },
            # Defaults
            defaults = {
                "general.default_runner": "SimpleJobRunner"
            },
            # Fallbacks
            fallbacks = {
                "runners.*": "general.default_runner"
            },
            # Configuration file
            settings_file=self._find_config(settings_file),
            resolve_undefined=resolve_undefined)

    def _find_config(self, settings_file=None):
        """
        Locate settings file to load values from

        Arguments:
          settings_file (str): if set then use this
            as the settings file; otherwise attempt
            locate a file by searching the default
            locations
        """
        if settings_file is None:
            # Look for default
            settings_file = locate_settings_file(
                "bcf_nanopore.ini",
                paths=[os.getcwd(),
                       get_config_dir()])
        return settings_file
