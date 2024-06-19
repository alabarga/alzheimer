# This snippet illustrates how to query David from R, using the RDAVIDWebService package.

# Load RDAVIDWebService.
library("RDAVIDWebService")

# Create a DAVIDWebService object connected to David, using your registration email.
# To register, go to: http://david.abcc.ncifcrf.gov/content.jsp?file=WS.html.
david <- DAVIDWebService$new(email='yourEmail@yourEmail.org')

# Define foreground and background gene lists.
# The foreground list should be contained within the background list.
FG <- addList(david, myForegroundGenes, idType="REFSEQ_MRNA", listName="isClass", listType="Gene")
BG <- addList(david, myBackgroundGenes, idType="REFSEQ_MRNA", listName="all", listType="Background")

# Inspect FG and BG to see the proportion of genes recognized by David, and those that are unmapped.
FG
BG

# Inspect "david" object to see the gene lists selected as foreground and background.
david

# Specifiy annotation categories.
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))

# Get functional annotation chart as R object.
FuncAnnotChart <- getFunctionalAnnotationChart(david)

# Print functional annotation chart to file.
getFunctionalAnnotationChartFile(david, "FuncAnnotChart.tsv")

# Get functional annotation clustering (limited to 3000 genes).
FuncAnnotClust <- getClusterReport(david)

# Print functional annotation clustering to file (limited to 3000 genes).
getClusterReportFile(david, "FuncAnnotClust.tsv")