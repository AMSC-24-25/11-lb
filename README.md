# Hands-on Group 11: Lattice Boltzmann Methods 
### Introduction and objectives
This projects uses the Lattice Boltzmann Methods (LBM) to perform a 2D fluid simulation with the D2Q9 model. Specifically, the aim was to solve the 2D lid-driven cavity problem.


### Building the Project

To build the project using CMake, follow these steps:

1. Clone the repository and navigate to the main directory.
2. Run the following commands:

   
   ### Create a directory for build files
   ```bash
   mkdir build
   ```

   
   ### Enter the build directory
   ```bash
   cd build
   ```
  
   ### Generate configuration files
   ```bash
   cmake ..
   ```
    
   ### Build the project
   ```bash
   cmake --build .
   ```



### Parallelization Strategy
We parallelized the code using OpenMP. The main bottleneck in the computation is represented by a single nested for-loop. The code was written so that each iteration was completely independent from the others and so that they all could be executed in parallel. To achieve this we needed to *bufferize* the whole computation so that two sets of memory location were used and then swapped for each iteration.
We observed a significant improvement in computation time.


### Scalability

<div class="container">
  <div class="column">
    <h4>Strong scalability test (6 cores, 12 threads)</h4>
    <div class="note">Execution times recorded on a PC with 6 cores and 12 threads.<br>Simulation run for 10'000 iterations.</div>
    <table>
      <thead>
        <tr>
          <th>Threads</th>
          <th>Time (ms)</th>
          <th>Speedup</th>
          <th>Efficiency</th>
        </tr>
      </thead>
      <tbody>
        <tr><td>1</td><td>69187</td><td>-</td><td>-</td></tr>
        <tr><td>2</td><td>38153</td><td>1.81</td><td>0.91</td></tr>
        <tr><td>3</td><td>29987</td><td>2.31</td><td>0.77</td></tr>
        <tr><td>4</td><td>27700</td><td>2.50</td><td>0.63</td></tr>
        <tr><td>5</td><td>24697</td><td>2.80</td><td>0.56</td></tr>
        <tr><td>6</td><td>22624</td><td>3.06</td><td>0.51</td></tr>
        <tr><td>7</td><td>19980</td><td>3.46</td><td>0.49</td></tr>
        <tr><td>8</td><td>17710</td><td>3.91</td><td>0.49</td></tr>
        <tr><td>9</td><td>16086</td><td>4.30</td><td>0.48</td></tr>
        <tr><td>10</td><td>14747</td><td>4.69</td><td>0.47</td></tr>
        <tr><td>11</td><td>14664</td><td>4.72</td><td>0.43</td></tr>
        <tr><td>12</td><td>18187</td><td>3.80</td><td>0.32</td></tr>
      </tbody>
    </table>
  </div>
  
  <div class="column">
    <h4>Strong scalability test (8 cores, 8 threads)</h4>
    <div class="note">Execution times recorded on a PC with 8 cores and 8 threads.<br>Simulation run for 10'000 iterations.</div>
    <table>
      <thead>
        <tr>
          <th>Threads</th>
          <th>Time (ms)</th>
          <th>Speedup</th>
          <th>Efficiency</th>
        </tr>
      </thead>
      <tbody>
        <tr><td>1</td><td>53929</td><td>-</td><td>-</td></tr>
        <tr><td>2</td><td>31645</td><td>1.70</td><td>0.85</td></tr>
        <tr><td>3</td><td>23348</td><td>2.31</td><td>0.77</td></tr>
        <tr><td>4</td><td>19208</td><td>2.81</td><td>0.70</td></tr>
        <tr><td>5</td><td>16473</td><td>3.27</td><td>0.65</td></tr>
        <tr><td>6</td><td>15139</td><td>3.56</td><td>0.59</td></tr>
        <tr><td>7</td><td>13518</td><td>3.99</td><td>0.57</td></tr>
        <tr><td>8</td><td>25401</td><td>2.12</td><td>0.27</td></tr>
      </tbody>
    </table>
  </div>
</div>


<div class="container">
  <div class="column">
    <h4>Weak scalability test (6 cores, 12 threads)</h4>
    <div class="note">Execution times recorded on a PC with 6 cores and 12 threads.<br>Simulation run for 1'000 iterations with ~20'000 lattice points per thread.</div>
    <table>
      <thead>
        <tr>
          <th>Threads</th>
          <th>Cavity Size</th>
          <th>Time (ms)</th>
          <th>Efficiency</th>
        </tr>
      </thead>
      <tbody>
        <tr><td>1</td><td>142x141</td><td>1381</td><td>-</td></tr>
        <tr><td>2</td><td>200x200</td><td>1537</td><td>0.90</td></tr>
        <tr><td>3</td><td>246x244</td><td>1821</td><td>0.76</td></tr>
        <tr><td>4</td><td>284x282</td><td>2029</td><td>0.68</td></tr>
        <tr><td>5</td><td>317x316</td><td>2205</td><td>0.63</td></tr>
        <tr><td>6</td><td>347x346</td><td>2351</td><td>0.59</td></tr>
        <tr><td>7</td><td>375x374</td><td>2503</td><td>0.55</td></tr>
        <tr><td>8</td><td>400x400</td><td>2598</td><td>0.53</td></tr>
        <tr><td>9</td><td>425x424</td><td>2651</td><td>0.52</td></tr>
        <tr><td>10</td><td>448x447</td><td>2748</td><td>0.50</td></tr>
        <tr><td>11</td><td>470x469</td><td>2942</td><td>0.47</td></tr>
        <tr><td>12</td><td>491x489</td><td>3288</td><td>0.42</td></tr>
      </tbody>
    </table>
  </div>
  
  <div class="column">
    <h4>Weak scalability test (8 cores, 8 threads)</h4>
    <div class="note">Execution times recorded on a PC with 8 cores and 8 threads.<br>Simulation run for 1'000 iterations with ~20'000 lattice points per thread.</div>
    <table>
      <thead>
        <tr>
          <th>Threads</th>
          <th>Cavity Size</th>
          <th>Time (ms)</th>
          <th>Efficiency</th>
        </tr>
      </thead>
      <tbody>
        <tr><td>1</td><td>142x141</td><td>1016</td><td>-</td></tr>
        <tr><td>2</td><td>200x200</td><td>1027</td><td>0.99</td></tr>
        <tr><td>3</td><td>246x244</td><td>1058</td><td>0.96</td></tr>
        <tr><td>4</td><td>284x282</td><td>1109</td><td>0.92</td></tr>
        <tr><td>5</td><td>317x316</td><td>1191</td><td>0.85</td></tr>
        <tr><td>6</td><td>347x346</td><td>1257</td><td>0.81</td></tr>
        <tr><td>7</td><td>375x374</td><td>1295</td><td>0.78</td></tr>
        <tr><td>8</td><td>400x400</td><td>1413</td><td>0.72</td></tr>
      </tbody>
    </table>
  </div>
</div>

### Validation of Results
We were provided with some reference data to compare our results.
The following tables, for instance, allows to compare the y-component of the velocity of the fluid along the vertical line through the geometrical center of the cavity.  

It can be observed that although the results with Re = 100 are relatively accurate, there's an error build-up with toward the lower part of the cavity and with increasing Reynolds numbers. That can probably be caused by slightly different parameters and formulae used. We are confident that our code works as intended and that, with some tweaking, we could allign our results with the reference data.



<div class="container2">
  <div class="row2">
    <div class="column2">
      <h4>Our results</h4>
      <table>
        <thead>
          <tr>
            <th>129-grid pt. no.</th>
            <th>y</th>
            <th>Re (100)</th>
            <th>Re (400)</th>
            <th>Re (1000)</th>
          </tr>
        </thead>
        <tbody>
          <tr><td>129</td><td>1.00</td><td>1</td><td>1</td><td>1</td></tr>
          <tr><td>126</td><td>0.98</td><td>0.844576</td><td>0.759886</td><td>0.654879</td></tr>
          <tr><td>125</td><td>0.97</td><td>0.792432</td><td>0.68533</td><td>0.56604</td></tr>
          <tr><td>124</td><td>0.96</td><td>0.741241</td><td>0.6171</td><td>0.496168</td></tr>
          <tr><td>123</td><td>0.95</td><td>0.691354</td><td>0.556289</td><td>0.44422</td></tr>
          <tr><td>110</td><td>0.85</td><td>0.222347</td><td>0.262017</td><td>0.295701</td></tr>
          <tr><td>95</td><td>0.74</td><td>-0.01649</td><td>0.141458</td><td>0.163419</td></tr>
          <tr><td>80</td><td>0.62</td><td>-0.1383</td><td>0.004383</td><td>0.042401</td></tr>
          <tr><td>65</td><td>0.50</td><td>-0.17906</td><td>-0.13328</td><td>-0.0651</td></tr>
          <tr><td>59</td><td>0.46</td><td>-0.17623</td><td>-0.18906</td><td>-0.10638</td></tr>
          <tr><td>37</td><td>0.29</td><td>-0.12114</td><td>-0.27501</td><td>-0.2716</td></tr>
          <tr><td>23</td><td>0.18</td><td>-0.07795</td><td>-0.16957</td><td>-0.31972</td></tr>
          <tr><td>14</td><td>0.11</td><td>-0.04994</td><td>-0.09319</td><td>-0.20474</td></tr>
          <tr><td>10</td><td>0.08</td><td>-0.03668</td><td>-0.06236</td><td>-0.14069</td></tr>
          <tr><td>9</td><td>0.07</td><td>-0.03321</td><td>-0.05479</td><td>-0.1247</td></tr>
          <tr><td>8</td><td>0.06</td><td>-0.02968</td><td>-0.04723</td><td>-0.10798</td></tr>
          <tr><td>1</td><td>0.00</td><td>0</td><td>0</td><td>0</td></tr>
        </tbody>
      </table>
    </div>
    <div class="column2">
      <h4>Reference data</h4>
      <table>
        <thead>
          <tr>
            <th>129-grid pt. no.</th>
            <th>y</th>
            <th>Re (100)</th>
            <th>Re (400)</th>
            <th>Re (1000)</th>
          </tr>
        </thead>
        <tbody>
          <tr><td>129</td><td>1.00000</td><td>1.00000</td><td>1.00000</td><td>1.00000</td></tr>
          <tr><td>126</td><td>0.9766</td><td>0.84123</td><td>0.75837</td><td>0.65928</td></tr>
          <tr><td>125</td><td>0.9688</td><td>0.78871</td><td>0.68439</td><td>0.57492</td></tr>
          <tr><td>124</td><td>0.9609</td><td>0.73722</td><td>0.61756</td><td>0.51117</td></tr>
          <tr><td>123</td><td>0.9531</td><td>0.68717</td><td>0.55892</td><td>0.46604</td></tr>
          <tr><td>110</td><td>0.8516</td><td>0.23151</td><td>0.29093</td><td>0.33304</td></tr>
          <tr><td>95</td><td>0.7344</td><td>0.00332</td><td>0.16256</td><td>0.18719</td></tr>
          <tr><td>80</td><td>0.6172</td><td>-0.13641</td><td>0.02135</td><td>0.05702</td></tr>
          <tr><td>65</td><td>0.5000</td><td>-0.20581</td><td>-0.11477</td><td>-0.06080</td></tr>
          <tr><td>59</td><td>0.4531</td><td>-0.21090</td><td>-0.17119</td><td>-0.10648</td></tr>
          <tr><td>37</td><td>0.2813</td><td>-0.15662</td><td>-0.32726</td><td>-0.27805</td></tr>
          <tr><td>23</td><td>0.1719</td><td>-0.10150</td><td>-0.24299</td><td>-0.38289</td></tr>
          <tr><td>14</td><td>0.1016</td><td>-0.06434</td><td>-0.14612</td><td>-0.29730</td></tr>
          <tr><td>10</td><td>0.0703</td><td>-0.04775</td><td>-0.10388</td><td>-0.22220</td></tr>
          <tr><td>9</td><td>0.0625</td><td>-0.04192</td><td>-0.09266</td><td>-0.20196</td></tr>
          <tr><td>8</td><td>0.0547</td><td>-0.03717</td><td>-0.08186</td><td>-0.18109</td></tr>
          <tr><td>1</td><td>0.0000</td><td>0.00000</td><td>0.00000</td><td>0.00000</td></tr>
        </tbody>
      </table>
    </div>
  </div>
  
  <div class="centered2">
    <div style="width: 70%;">
      <h4>Error wrt the reference data</h4>
      <table>
        <thead>
          <tr>
            <th>129-Grid pt. no.</th>
            <th>y</th>
            <th>Re (100)</th>
            <th>Re (400)</th>
            <th>Re (1000)</th>
          </tr>
        </thead>
        <tbody>
          <tr><td>129</td><td>1.00</td><td>0.000000</td><td>0.000000</td><td>0.000000</td></tr>
          <tr><td>126</td><td>0.98</td><td>0.003346</td><td>0.001516</td><td>-0.004401</td></tr>
          <tr><td>125</td><td>0.97</td><td>0.003722</td><td>0.000940</td><td>-0.008880</td></tr>
          <tr><td>124</td><td>0.96</td><td>0.004021</td><td>-0.000460</td><td>-0.015002</td></tr>
          <tr><td>123</td><td>0.95</td><td>0.004184</td><td>-0.002631</td><td>-0.021820</td></tr>
          <tr><td>110</td><td>0.85</td><td>-0.009163</td><td>-0.028913</td><td>-0.037339</td></tr>
          <tr><td>95</td><td>0.74</td><td>-0.019810</td><td>-0.021102</td><td>-0.023771</td></tr>
          <tr><td>80</td><td>0.62</td><td>-0.001890</td><td>-0.016967</td><td>-0.014619</td></tr>
          <tr><td>65</td><td>0.50</td><td>0.026750</td><td>-0.018510</td><td>-0.004300</td></tr>
          <tr><td>59</td><td>0.46</td><td>0.034670</td><td>-0.017870</td><td>0.000100</td></tr>
          <tr><td>37</td><td>0.29</td><td>0.035480</td><td>0.052250</td><td>0.006450</td></tr>
          <tr><td>23</td><td>0.18</td><td>0.023550</td><td>0.073420</td><td>0.063170</td></tr>
          <tr><td>14</td><td>0.11</td><td>0.014400</td><td>0.052930</td><td>0.092560</td></tr>
          <tr><td>10</td><td>0.08</td><td>0.011070</td><td>0.041520</td><td>0.081510</td></tr>
          <tr><td>9</td><td>0.07</td><td>0.008710</td><td>0.037870</td><td>0.077260</td></tr>
          <tr><td>8</td><td>0.06</td><td>0.007490</td><td>0.034630</td><td>0.073110</td></tr>
          <tr><td>1</td><td>0.00</td><td>0.000000</td><td>0.000000</td><td>0.000000</td></tr>
        </tbody>
      </table>
    </div>
  </div>
</div>


### Gallery

<table>
    <tr>
        <th>Re=100 on a 100x100 cavity</th>
        <th>Re=10000 on a 100x100 cavity</th>
    </tr>
    <tr>
        <td><img src="./media/100x100_re100_steps5000_periteration5_fps24.gif" alt="100x100_re100_steps5000_periteration5_fps24"></td>
        <td><img src="./media/100x100_re10000_steps10000_periteration10_fps24.gif" alt="100x100_re10000_steps10000_periteration10_fps24"></td>
    </tr>
    <tr>
        <th>Re=10000 on a 75x50 cavity with an alternative boundary condition</th>
        <th>Re=10000 on a 240x200 cavity with an alternative boundary condition</th>
    </tr>
    <tr>
        <td><img src="./media/75x50_re10000_steps10000_periteration20_fps24_alt1.gif" alt="75x50_re10000_steps10000_periteration20_fps24_alt1"></td>
        <td><img src="./media/240x200_re10000_steps20000_periteration50_fps24_alt2.gif" alt="240x200_re10000_steps20000_periteration50_fps24_alt2"></td>
    </tr>
    <tr>
        <th>Re=2200 on a 200x80 cavity </th>
        <th>Re=1000 on a 250x250 cavity</th>
    </tr>
    <tr>
        <td><img src="./media/200x80_re2200_steps10000_periteration50_fps24.gif" alt="200x80_re2200_steps10000_periteration50_fps24"></td>
        <td><img src="./media/250x250_re1000_steps10000_periteration20_fps24.gif" alt="250x250_re1000_steps10000_periteration20_fps24"></td>
    </tr>
</table>


