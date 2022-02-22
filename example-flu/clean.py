#!/usr/bin/env python3

import os
import ppdx

import logging
log = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(process)s - %(message)s')

def main():

    ppdx.config.cread('config-ppdx.ini')
    ppdx.WRKDIR = os.path.join(os.getcwd(), "models")
    ppdx.clean()
    log.info("Finished!")

main()

