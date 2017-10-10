# TReNA
 Fit transcriptional regulatory networks using gene expression, priors, machine learning

## Getting Started 

To build and test:

 - clone this repository
 - install R 3.2.3 or later; RUnit 0.4.31 or later (see below)
 - install the following solver packages:
   - glmnet R package 2.0.3 or later
   - randomForest
   - vbsr
   - flare
   - lassopv
 - cd TReNA
 - R CMD INSTALL .
 
The most reliable way to install package dependencies (and other of their dependencies):

````
source("http://bioconductor.org/biocLite.R")
biocLite(c("glmnet", "RUnit"))
````

## Using TReNA

 - open an R session
 - source("inst/unitTests/test_TReNA.R")
 - runTests()

The unitTests perform double duty: they ensure the package performs as (currently) expected;
they introduce the package to the user and developer.
Thus [test_TReNA.R](https://github.com/PriceLab/TReNA/blob/master/inst/unitTests/test_TReNA.R)
is one entry point into this project.

We have also created a [Jupyter Notebook](http://nbviewer.jupyter.org/github/PriceLab/TReNA/blob/master/inst/demos/Assess_Distributions.ipynb) demonstrating use of TReNA with 4 different solvers

## Table of BDDS Databases (2017/10/10)

<table class='gmisc_table' style='border-collapse: collapse; margin-top: 1em; margin-bottom: 1em;' >
<thead>
<tr>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey;'> </th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Name</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Hits</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Hits.thousand</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Hits.million</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Tissue</th>
<th style='border-bottom: 1px solid grey; border-top: 2px solid grey; text-align: center;'>Method</th>
</tr>
</thead>
<tbody>
<tr>
<td style='text-align: left;'>1</td>
<td style='text-align: center;'>adrenal_gland_hint_20</td>
<td style='text-align: center;'>33500100</td>
<td style='text-align: center;'>33500.1</td>
<td style='text-align: center;'>33.5001</td>
<td style='text-align: center;'>adrenal_gland</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>2</td>
<td style='text-align: center;'>adrenal_gland_hint_16</td>
<td style='text-align: center;'>81369500</td>
<td style='text-align: center;'>81369.5</td>
<td style='text-align: center;'>81.3695</td>
<td style='text-align: center;'>adrenal_gland</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>3</td>
<td style='text-align: center;'>bone_element_hint_20</td>
<td style='text-align: center;'>33500100</td>
<td style='text-align: center;'>33500.1</td>
<td style='text-align: center;'>33.5001</td>
<td style='text-align: center;'>bone_element</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>4</td>
<td style='text-align: center;'>brain_hint_16</td>
<td style='text-align: center;'>839776000</td>
<td style='text-align: center;'>839776</td>
<td style='text-align: center;'>839.776</td>
<td style='text-align: center;'>brain</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>5</td>
<td style='text-align: center;'>brain_hint_20</td>
<td style='text-align: center;'>865315000</td>
<td style='text-align: center;'>865315</td>
<td style='text-align: center;'>865.315</td>
<td style='text-align: center;'>brain</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>6</td>
<td style='text-align: center;'>bronchus_hint_20</td>
<td style='text-align: center;'>33500100</td>
<td style='text-align: center;'>33500.1</td>
<td style='text-align: center;'>33.5001</td>
<td style='text-align: center;'>bronchus</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>7</td>
<td style='text-align: center;'>bronchus_hint_16</td>
<td style='text-align: center;'>62819300</td>
<td style='text-align: center;'>62819.3</td>
<td style='text-align: center;'>62.8193</td>
<td style='text-align: center;'>bronchus</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>8</td>
<td style='text-align: center;'>esophagus_hint_20</td>
<td style='text-align: center;'>55235500</td>
<td style='text-align: center;'>55235.5</td>
<td style='text-align: center;'>55.2355</td>
<td style='text-align: center;'>esophagus</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>9</td>
<td style='text-align: center;'>extraembryonic_structure_hint_16</td>
<td style='text-align: center;'>386228000</td>
<td style='text-align: center;'>386228</td>
<td style='text-align: center;'>386.228</td>
<td style='text-align: center;'>extraembryonic_structure</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>10</td>
<td style='text-align: center;'>extraembryonic_structure_hint_20</td>
<td style='text-align: center;'>391720000</td>
<td style='text-align: center;'>391720</td>
<td style='text-align: center;'>391.72</td>
<td style='text-align: center;'>extraembryonic_structure</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>11</td>
<td style='text-align: center;'>eye_hint_20</td>
<td style='text-align: center;'>322568000</td>
<td style='text-align: center;'>322568</td>
<td style='text-align: center;'>322.568</td>
<td style='text-align: center;'>eye</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>12</td>
<td style='text-align: center;'>eye_hint_16</td>
<td style='text-align: center;'>324257000</td>
<td style='text-align: center;'>324257</td>
<td style='text-align: center;'>324.257</td>
<td style='text-align: center;'>eye</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>13</td>
<td style='text-align: center;'>gonad_hint_20</td>
<td style='text-align: center;'>39931800</td>
<td style='text-align: center;'>39931.8</td>
<td style='text-align: center;'>39.9318</td>
<td style='text-align: center;'>gonad</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>14</td>
<td style='text-align: center;'>heart_hint_20</td>
<td style='text-align: center;'>295442000</td>
<td style='text-align: center;'>295442</td>
<td style='text-align: center;'>295.442</td>
<td style='text-align: center;'>heart</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>15</td>
<td style='text-align: center;'>kidney_hint_20</td>
<td style='text-align: center;'>289504000</td>
<td style='text-align: center;'>289504</td>
<td style='text-align: center;'>289.504</td>
<td style='text-align: center;'>kidney</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>16</td>
<td style='text-align: center;'>large_intestine_hint_20</td>
<td style='text-align: center;'>125992000</td>
<td style='text-align: center;'>125992</td>
<td style='text-align: center;'>125.992</td>
<td style='text-align: center;'>large_intestine</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>17</td>
<td style='text-align: center;'>liver_hint_20</td>
<td style='text-align: center;'>63053400</td>
<td style='text-align: center;'>63053.4</td>
<td style='text-align: center;'>63.0534</td>
<td style='text-align: center;'>liver</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>18</td>
<td style='text-align: center;'>lung_hint_20</td>
<td style='text-align: center;'>252351000</td>
<td style='text-align: center;'>252351</td>
<td style='text-align: center;'>252.351</td>
<td style='text-align: center;'>lung</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>19</td>
<td style='text-align: center;'>lymphatic_vessel_hint_20</td>
<td style='text-align: center;'>58418900</td>
<td style='text-align: center;'>58418.9</td>
<td style='text-align: center;'>58.4189</td>
<td style='text-align: center;'>lymphatic_vessel</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>20</td>
<td style='text-align: center;'>lymphoblast_hint_16</td>
<td style='text-align: center;'>483616000</td>
<td style='text-align: center;'>483616</td>
<td style='text-align: center;'>483.616</td>
<td style='text-align: center;'>lymphoblast</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>21</td>
<td style='text-align: center;'>lymphoblast_hint_20</td>
<td style='text-align: center;'>498033000</td>
<td style='text-align: center;'>498033</td>
<td style='text-align: center;'>498.033</td>
<td style='text-align: center;'>lymphoblast</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>22</td>
<td style='text-align: center;'>mammary_gland_hint_20</td>
<td style='text-align: center;'>58580200</td>
<td style='text-align: center;'>58580.2</td>
<td style='text-align: center;'>58.5802</td>
<td style='text-align: center;'>mammary_gland</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>23</td>
<td style='text-align: center;'>mouth_hint_20</td>
<td style='text-align: center;'>156634000</td>
<td style='text-align: center;'>156634</td>
<td style='text-align: center;'>156.634</td>
<td style='text-align: center;'>mouth</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>24</td>
<td style='text-align: center;'>muscle_organ_hint_20</td>
<td style='text-align: center;'>116657000</td>
<td style='text-align: center;'>116657</td>
<td style='text-align: center;'>116.657</td>
<td style='text-align: center;'>muscle_organ</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>25</td>
<td style='text-align: center;'>pancreas_hint_20</td>
<td style='text-align: center;'>81055800</td>
<td style='text-align: center;'>81055.8</td>
<td style='text-align: center;'>81.0558</td>
<td style='text-align: center;'>pancreas</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>26</td>
<td style='text-align: center;'>prostate_gland_hint_20</td>
<td style='text-align: center;'>49664100</td>
<td style='text-align: center;'>49664.1</td>
<td style='text-align: center;'>49.6641</td>
<td style='text-align: center;'>prostate_gland</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>27</td>
<td style='text-align: center;'>skin_hint_16</td>
<td style='text-align: center;'>1281170000</td>
<td style='text-align: center;'>1281170</td>
<td style='text-align: center;'>1281.17</td>
<td style='text-align: center;'>skin</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>28</td>
<td style='text-align: center;'>skin_hint_20</td>
<td style='text-align: center;'>1317110000</td>
<td style='text-align: center;'>1317110</td>
<td style='text-align: center;'>1317.11</td>
<td style='text-align: center;'>skin</td>
<td style='text-align: center;'>hint</td>
</tr>
<tr>
<td style='text-align: left;'>29</td>
<td style='text-align: center;'>adrenal_gland_wellington_16</td>
<td style='text-align: center;'>24344500</td>
<td style='text-align: center;'>24344.5</td>
<td style='text-align: center;'>24.3445</td>
<td style='text-align: center;'>adrenal_gland</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>30</td>
<td style='text-align: center;'>adrenal_gland_wellington_20</td>
<td style='text-align: center;'>24769200</td>
<td style='text-align: center;'>24769.2</td>
<td style='text-align: center;'>24.7692</td>
<td style='text-align: center;'>adrenal_gland</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>31</td>
<td style='text-align: center;'>bone_element_wellington_20</td>
<td style='text-align: center;'>7174280</td>
<td style='text-align: center;'>7174.28</td>
<td style='text-align: center;'>7.17428</td>
<td style='text-align: center;'>bone_element</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>32</td>
<td style='text-align: center;'>brain_wellington_16</td>
<td style='text-align: center;'>228432000</td>
<td style='text-align: center;'>228432</td>
<td style='text-align: center;'>228.432</td>
<td style='text-align: center;'>brain</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>33</td>
<td style='text-align: center;'>brain_wellington_20</td>
<td style='text-align: center;'>248069000</td>
<td style='text-align: center;'>248069</td>
<td style='text-align: center;'>248.069</td>
<td style='text-align: center;'>brain</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>34</td>
<td style='text-align: center;'>bronchus_wellington_20</td>
<td style='text-align: center;'>7174280</td>
<td style='text-align: center;'>7174.28</td>
<td style='text-align: center;'>7.17428</td>
<td style='text-align: center;'>bronchus</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>35</td>
<td style='text-align: center;'>bronchus_wellington_16</td>
<td style='text-align: center;'>13069200</td>
<td style='text-align: center;'>13069.2</td>
<td style='text-align: center;'>13.0692</td>
<td style='text-align: center;'>bronchus</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>36</td>
<td style='text-align: center;'>esophagus_wellington_20</td>
<td style='text-align: center;'>17559700</td>
<td style='text-align: center;'>17559.7</td>
<td style='text-align: center;'>17.5597</td>
<td style='text-align: center;'>esophagus</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>37</td>
<td style='text-align: center;'>extraembryonic_structure_wellington_16</td>
<td style='text-align: center;'>25254000</td>
<td style='text-align: center;'>25254</td>
<td style='text-align: center;'>25.254</td>
<td style='text-align: center;'>extraembryonic_structure</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>38</td>
<td style='text-align: center;'>extraembryonic_structure_wellington_20</td>
<td style='text-align: center;'>116007000</td>
<td style='text-align: center;'>116007</td>
<td style='text-align: center;'>116.007</td>
<td style='text-align: center;'>extraembryonic_structure</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>39</td>
<td style='text-align: center;'>eye_wellington_20</td>
<td style='text-align: center;'>81794200</td>
<td style='text-align: center;'>81794.2</td>
<td style='text-align: center;'>81.7942</td>
<td style='text-align: center;'>eye</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>40</td>
<td style='text-align: center;'>eye_wellington_16</td>
<td style='text-align: center;'>82913900</td>
<td style='text-align: center;'>82913.9</td>
<td style='text-align: center;'>82.9139</td>
<td style='text-align: center;'>eye</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>41</td>
<td style='text-align: center;'>gonad_wellington_20</td>
<td style='text-align: center;'>22184900</td>
<td style='text-align: center;'>22184.9</td>
<td style='text-align: center;'>22.1849</td>
<td style='text-align: center;'>gonad</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>42</td>
<td style='text-align: center;'>heart_wellington_20</td>
<td style='text-align: center;'>80855200</td>
<td style='text-align: center;'>80855.2</td>
<td style='text-align: center;'>80.8552</td>
<td style='text-align: center;'>heart</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>43</td>
<td style='text-align: center;'>kidney_wellington_20</td>
<td style='text-align: center;'>74294800</td>
<td style='text-align: center;'>74294.8</td>
<td style='text-align: center;'>74.2948</td>
<td style='text-align: center;'>kidney</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>44</td>
<td style='text-align: center;'>large_intestine_wellington_20</td>
<td style='text-align: center;'>55642800</td>
<td style='text-align: center;'>55642.8</td>
<td style='text-align: center;'>55.6428</td>
<td style='text-align: center;'>large_intestine</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>45</td>
<td style='text-align: center;'>liver_wellington_20</td>
<td style='text-align: center;'>21973800</td>
<td style='text-align: center;'>21973.8</td>
<td style='text-align: center;'>21.9738</td>
<td style='text-align: center;'>liver</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>46</td>
<td style='text-align: center;'>lung_wellington_20</td>
<td style='text-align: center;'>56271600</td>
<td style='text-align: center;'>56271.6</td>
<td style='text-align: center;'>56.2716</td>
<td style='text-align: center;'>lung</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>47</td>
<td style='text-align: center;'>lymphatic_vessel_wellington_20</td>
<td style='text-align: center;'>10199500</td>
<td style='text-align: center;'>10199.5</td>
<td style='text-align: center;'>10.1995</td>
<td style='text-align: center;'>lymphatic_vessel</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>48</td>
<td style='text-align: center;'>lymphoblast_wellington_16</td>
<td style='text-align: center;'>150023000</td>
<td style='text-align: center;'>150023</td>
<td style='text-align: center;'>150.023</td>
<td style='text-align: center;'>lymphoblast</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>49</td>
<td style='text-align: center;'>lymphoblast_wellington_20</td>
<td style='text-align: center;'>155832000</td>
<td style='text-align: center;'>155832</td>
<td style='text-align: center;'>155.832</td>
<td style='text-align: center;'>lymphoblast</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>50</td>
<td style='text-align: center;'>mammary_gland_wellington_20</td>
<td style='text-align: center;'>10972000</td>
<td style='text-align: center;'>10972</td>
<td style='text-align: center;'>10.972</td>
<td style='text-align: center;'>mammary_gland</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>51</td>
<td style='text-align: center;'>mouth_wellington_20</td>
<td style='text-align: center;'>40040100</td>
<td style='text-align: center;'>40040.1</td>
<td style='text-align: center;'>40.0401</td>
<td style='text-align: center;'>mouth</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>52</td>
<td style='text-align: center;'>muscle_organ_wellington_20</td>
<td style='text-align: center;'>30304400</td>
<td style='text-align: center;'>30304.4</td>
<td style='text-align: center;'>30.3044</td>
<td style='text-align: center;'>muscle_organ</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>53</td>
<td style='text-align: center;'>pancreas_wellington_20</td>
<td style='text-align: center;'>36360100</td>
<td style='text-align: center;'>36360.1</td>
<td style='text-align: center;'>36.3601</td>
<td style='text-align: center;'>pancreas</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>54</td>
<td style='text-align: center;'>prostate_gland_wellington_20</td>
<td style='text-align: center;'>17724100</td>
<td style='text-align: center;'>17724.1</td>
<td style='text-align: center;'>17.7241</td>
<td style='text-align: center;'>prostate_gland</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='text-align: left;'>55</td>
<td style='text-align: center;'>skin_wellington_20</td>
<td style='text-align: center;'>0</td>
<td style='text-align: center;'>0</td>
<td style='text-align: center;'>0</td>
<td style='text-align: center;'>skin</td>
<td style='text-align: center;'>wellington</td>
</tr>
<tr>
<td style='border-bottom: 2px solid grey; text-align: left;'>56</td>
<td style='border-bottom: 2px solid grey; text-align: center;'>skin_wellington_16</td>
<td style='border-bottom: 2px solid grey; text-align: center;'>287181000</td>
<td style='border-bottom: 2px solid grey; text-align: center;'>287181</td>
<td style='border-bottom: 2px solid grey; text-align: center;'>287.181</td>
<td style='border-bottom: 2px solid grey; text-align: center;'>skin</td>
<td style='border-bottom: 2px solid grey; text-align: center;'>wellington</td>
</tr>
</tbody>
</table>

