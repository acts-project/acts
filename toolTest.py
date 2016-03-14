###############################################################
#
# Job options
#
#==============================================================

from Gaudi.Configuration import *
from Configurables import ApplicationMgr

from ToolTest.ToolTestConf import ToolTest2
Tool2 = ToolTest2("Tool2")

from ToolTest.ToolTestConf import ToolTest
Tool1 = ToolTest("Tool1")
Tool1.ToolTest2 = Tool2

from ToolAlgorithm.ToolAlgorithmConf import ToolAlgorithm
ToolAlgorithm = ToolAlgorithm("ToolAlgorithm")
ToolAlgorithm.ToolTest = Tool1

ApplicationMgr(EvtSel='NONE',
               EvtMax=1,
               ExtSvc=[],
               TopAlg=[ToolAlgorithm])


#==============================================================
#
# End of job options file
#
###############################################################

