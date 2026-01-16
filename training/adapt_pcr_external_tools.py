#!/usr/bin/env python
# coding: utf-8
"""External tool path configuration for training scripts.

**Packages**
- RNAstructure: version 6.5 (released on Jun 14, 2024)
- primer3: libprimer3 release 2.6.1

Tool paths can be configured via environment variables:
- QPRIMER_TOOLPATH: Base path for tools (if tools are in a custom location)
- RNASTRUCTURE_DATAPATH: Path to RNAstructure data tables

If QPRIMER_TOOLPATH is not set, tools are assumed to be in PATH
(e.g., installed via conda).
"""

import os
import shutil


def _find_tool(name: str, toolpath: str | None = None, subpath: str | None = None) -> str:
    """Find a tool either in TOOLPATH or in system PATH."""
    if toolpath and subpath:
        full_path = os.path.join(toolpath, subpath)
        if os.path.isfile(full_path) and os.access(full_path, os.X_OK):
            return full_path

    # Fall back to PATH
    path = shutil.which(name)
    if path:
        return path

    raise FileNotFoundError(
        f"{name} not found. Install via conda or set QPRIMER_TOOLPATH environment variable."
    )


# Get optional tool path from environment
TOOLPATH = os.environ.get('QPRIMER_TOOLPATH')
DATAPATH = os.environ.get(
    'RNASTRUCTURE_DATAPATH',
    os.path.join(TOOLPATH, 'RNAstructure/data_tables') if TOOLPATH else ''
)

# Tool paths - use TOOLPATH if set, otherwise find in PATH
Fold = _find_tool('Fold', TOOLPATH, 'RNAstructure/exe/Fold')
bifold = _find_tool('bifold', TOOLPATH, 'RNAstructure/exe/bifold')
efn2 = _find_tool('efn2', TOOLPATH, 'RNAstructure/exe/efn2')
oligotm = _find_tool('oligotm', TOOLPATH, 'primer3/src/oligotm')
primer3 = _find_tool('primer3_core', TOOLPATH, 'primer3/src/primer3_core')
RNAduplex = _find_tool('RNAduplex', TOOLPATH, 'ViennaRNA-2.7.0/src/bin/RNAduplex')
bowtie2 = _find_tool('bowtie2', TOOLPATH, 'bowtie2-2.5.3-linux-x86_64/bowtie2')
sam2pairwise = _find_tool('sam2pairwise', TOOLPATH, 'sam2pairwise/src/sam2pairwise')


def show_available_tools():
    """Print available tools and their paths."""
    tools = {
        'Fold': Fold,
        'bifold': bifold,
        'efn2': efn2,
        'oligotm': oligotm,
        'primer3': primer3,
        'RNAduplex': RNAduplex,
        'bowtie2': bowtie2,
        'sam2pairwise': sam2pairwise,
    }
    print('Available tools:')
    for name, path in tools.items():
        print(f'  {name}: {path}')




