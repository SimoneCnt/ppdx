#!/usr/bin/env python3

import os
import ppdg

import logging
log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(process)s - %(message)s')

def main():

    ppdg.config.cread('config-ppdg.ini')
    ppdg.WRKDIR = os.path.join(os.getcwd(), "models")
    ppdg.clean()
    log.info("Finished!")

main()

