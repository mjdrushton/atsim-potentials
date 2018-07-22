import logging

from ...config import Configuration

def action_tabulate(cfg_file, outfilename):
  logger = logging.getLogger(__name__).getChild('_action_tabulate')
  config = Configuration()
  tabulation = config.read(cfg_file)
  with open(outfilename, 'w') as outfile:
    logger.info("Writing output to: {}".format(outfilename))
    tabulation.write(outfile)
