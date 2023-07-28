############################################################
# For personal site, 03.20
# This generates the file "index.html"
############################################################

import sys, os
sys.path.append(os.path.abspath('../lib/'))
import read_chunks as RC

######################
# HTML template
######################

html_template = """
<!doctype html>
    {head}

<body>
	<!-- {nav} 

    <div class="row" id="top_grid">
		<div class="col-2-24" id="margin"></div>
		<div class="col-20-24" id="main_header">Mus t-haplotype phylogenetics</div>
		<div class="col-2-24" id="margin"></div>
	</div>
    
	<div class="sep_div"></div>

	<div class="row">
		<div class="col-2-24" id="margin"></div>
		<div class="col-20-24" id="main_col">
			<ul>
				<li><h3><a href="https://github.com/gwct/mus-t-haplotype">github repo</a></h3></li>
				<li><h3><a href="cactus.html">Cactus stats</a><span style="font-size:0.65em;"></h3></li>
				<li><h3><a href="windows.html">Window stats</a><span style="font-size:0.65em;"></h3></li>
				<li><h3><a href="trees.html">Tree stats</a><span style="font-size:0.65em;"></h3></li>
			</ul>
		</div>
        <div class="col-2-24" id="margin"></div>
	</div>

	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
    <div class="sep_div"></div>
    <div class="sep_div"></div>

    {footer}
</body>
"""

######################
# Main block
######################
pagefile = "index.html";
print("Generating " + pagefile + "...");
title = "Mus t-haplotype"

head = RC.readHead(title, pagefile);
nav = RC.readNav(pagefile);
footer = RC.readFooter();

outfilename = "../../" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, footer=footer));