import logging

from ...config import Configuration

def action_tabulate(cp, outfilename):
  logger = logging.getLogger(__name__).getChild('_action_tabulate')
  config = Configuration()
  tabulation = config.read_from_parser(cp)
  with tabulation.open_fp(outfilename) as outfile:
    logger.info("Writing output to: {}".format(outfilename))
    tabulation.write(outfile)
