<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">

<head>
  <meta name="DC.Creator" content="Christiane Peters"></meta>
  <meta name="DC.Type" content="Text"></meta>
  <meta http-equiv="content-type" content="text/html; charset=utf-8"></meta>
  <meta http-equiv="Content-Style-Type" content="text/css"></meta>
  <link rel="stylesheet" type="text/css" href="style.css" />
  <title> Christiane Peters </title>
</head>

<body>
<h4>
Iteration and operation count for information-set decoding over F<sub>q</sub>
</h4>

<p> My paper <a href="./#2010.isdfq">Information-set decoding for linear codes over F<sub>q</sub></a> presents a new algorithm for decoding linear codes over arbitrary finite fields F<sub>q</sub>.
<br></br>
<br></br>
<b>Note</b>: if the field size q is large then the cost of the first step
(updating the matrix) becomes negligible since the cost of the algorithm is
dominated by the search for 2p columns adding up to 0
on l positions.
One can afford to choose k columns uniformly at random at the beginning of each iteration and then perform a full Gaussian elimination with respect to those 
k columns.

In order to get a good estimate for the running time you can use
<a href="isdfq.gp">this PARI/GP script</a> to estimate the cost of information-set decoding (ISD). <br></br>
Note that this rather crude approximation gives a running time which is worse
by a factor of 2 or more than the Markov-chain analysis below. 
It's fine if you just want to know in which ballpark your parameters are.

<br></br>
<br></br>
If q is small (q=2, q=3, maybe some cases where q=4)
it makes sense to speed up the first step (definitely for q=2, see  
<a href="./#2008.mceliece">our PQCrypto'08 paper</a>).
<br></br>
<br></br>
The ISD algorithm is probabilistic in that it makes random choices. 
To compute the average number of iterations it is necessary to analyze a Markov chain. 
For each choice of parameters of the decoding algorithm the following 
C program computes the expected number of iterations as well as the 
number of bit operations needed. 

If this approach is used to estimate the security of a code-based
cryptosystem against information-set-decoding attacks, the user must 
run the program  repeatedly to search a range of parameter choices 
for the decoding algorithm. For example, trying this program with a wide
range of algorithm parameters for decoding 48 errors in a [961,771] code 
over F<sub>31</sub> identifies the parameters p=2, l=7, c=12, and r=1 as optimal
and shows that 2<sup>128</sup> bit operations are used by the 
information-set-decoding attack, as reported in the paper.

In fact, it was used to search a range of attack parameters against a
range of codes over F<sub>31</sub> to identify this code with [961,771,>=96] as 
having minimal keysize k x (n-k) among all codes over F<sub>31</sub> that require 
at least 2<sup>128</sup> steps in the decoding algorithm.

See below for an explanation of how to use the program and several
examples.
</p>

<p>
<a href="isdfq.c">http://www2.mat.dtu.dk/people/C.Peters/isdfq.c</a>
</p>

<ul>
  <li>This program uses the
  <a href="http://perso.ens-lyon.fr/nathalie.revol/mpfi.html">
  MPFI library</a> which is built on top of the 
  <a href="http://www.mpfr.org/">MPFR library</a>,
  which is built on top of the 
  <a href="http://gmplib.org/">GMP library</a>.
  </li>
  <li>
  This program works for any field F<sub>q</sub>, including the binary case.
  This program supersedes the program given for binary codes presented at
  <a href="http://www2.mat.dtu.dk/people/C.Peters/mceliece.html">http://www2.mat.dtu.dk/people/C.Peters/mceliece.html</a>. 
  </li>
  <li>This program takes the following 13 inputs:
    <ul class="disc">
    <li>field size <tt> q </tt> (default 31)</li>
    <li>code length <tt> n </tt> (default 961); </li>
    <li>code dimension <tt> k </tt> (default 771); </li>
    <li>number of errors <tt> w </tt> (default 48); </li>
    <li>algorithm parameter <tt> p </tt> (default 2); </li>
    <li>algorithm parameter <tt> l </tt> (default 20); </li>
    <li>algorithm parameter <tt> m </tt> (default 1);  &nbsp; <a
	href="#m">(*)</a></li>
    <li>algorithm parameter <tt> c </tt> (default 7); </li>
    <li>algorithm parameter <tt> r </tt> (default 1); </li>
    <li>use overlapping sets by setting <tt> fs=1; </tt> (default 0); </li>
    <li>algorithm parameter <tt> M </tt> for the <tt>fs=1</tt> case (adjust number of subsets by multiplying the standard choice binomial(k,p)/sqrt(binomial(2p,p)) by a factor of M) (default 1.); </li>
    <li>look for a weight-w word by setting <tt> mww=1; </tt> (default 0); </li>
    <li>adjust precision by setting <tt> prec </tt> (default 300). </li>
    </ul>
  </li>
  <li>The user can confirm a choice of parameters for the algorithm using this program. The program does not choose the best algorithm parameters itself. The user should loop over all suitable choices for p,l,m,c,r, and M in order to figure out the lowest cost for an attack.
  </li>
  <li>Examples: 
    <ul class="disc">
    <li> Compile using <tt> gcc -o isdfq isdfq.c -lm -lgmp -lmpfr -lmpfi </tt></li>
    <li> Run <tt> ./isdfq </tt> to get the number of operations for decoding 48 errors in a [961,771] code over F<sub>31</sub> with parameters p=2, l=7, c=12, r=1:<br></br> 
    <ul class="none">
      <li> <tt>
        q=31 n=961 k=771 w=48 p=2 l=7 m=1 c=12 r=1:  
        bit ops  129.023892, bit ops per it 32.208747, 
        log2 #it 96.815146</tt>
      </li>
      </ul>
    </li>
    <li> Run <tt>./isdfq 4 2560 1536 128 2 15 1 8 4 0 1 0 400</tt>
    to get the number of operations it takes to find a weight-128 word in a 
    [2560,1536] code over F4 when using parameters
    p=2, l=15, m=1, c=8, and r=4.
    <ul class="none">
      <li><tt>
        q=4 n=2560 k=1536 w=128 p=2 l=15 m=1 c=8 r=4:  
        bit ops  181.857556, bit ops per it 27.515332, 
        log2 #it 154.342224</tt>
      </li>
    </ul>
    </li>
    <li> Run <tt> ./isdfq 2 1024 524 50 2 20 2 7 7</tt> to get the numbers of 
      operations to decode 50 errors in a [1024,524] binary code using parameters p=2, l=20, m=2, c=7, r=7.
      <ul class="none">
      <li><tt>q=2 n=1024 k=524 w=50 p=2 l=20 m=2 c=7 r=7:  
              bit ops  60.179978, bit ops per it 21.864968, 
              log2 #it 38.315010</tt>
      </li>
      </ul>
	</li>
    <li> Run <tt>./isdfq 2 1024 525 50 2 20 2 7 7 0 1. 1</tt>
         to get the number of operations it takes to find a weight-50 word in a 
         [1024,525] binary code using the algorithm with parameters 
         p=2, l=20, m=2, c=7, r=7 (note that the value of <tt>mww</tt> is set to 1).
    <ul class="none">
      <li><tt>q=2 n=1024 k=525 w=50 p=2 l=20 m=2 c=7 r=7:  
              bit ops  60.392514, bit ops per it 21.944141, 
              log2 #it 38.448373
          </tt>
      </li>
    </ul>
    </li>
  </ul>
  </li>
  <li> Cost of attacking codes proposed by 
    <ul>
      <li>Thierry Berger, Pierre-Louis Cayrel, Philippe Gaborit, Ayoub Otmani.
          <i>Reducing key length of the McEliece cryptosystem</i>. Africacrypt 2009.
	   <a href="http://users.info.unicaen.fr/~otmani/pdf/africacrypt09.pdf">
http://users.info.unicaen.fr/~otmani/pdf/africacrypt09.pdf</a>
      </li>
      <li>Rafael Misoczki, Paulo S. L. M. Barreto. <i>Compact McEliece keys from Goppa codes</i>. SAC 2009. 
       <a href="http://eprint.iacr.org/2009/187">http://eprint.iacr.org/2009/187</a>

      </li>
    </ul>

   <table border="1">

    <tr>
      <th colspan="4">code parameters</th>
      <th rowspan="2">claimed <br></br>security level</th>
      <th colspan="6">disjoint split</th>
      <th colspan="8">overlapping sets</th>
    </tr>
    <tr>
      <th>q </th>
      <th>n </th>
      <th>k </th>
      <th>w </th>
      <th>log2(# bit ops)</th>
   <!--   <th>log2(# bit ops/it)</th>
      <th>log2(# iterations)</th>-->
      <th>p </th>
      <th>l </th>
      <th>m </th>
      <th>c </th>
      <th>r </th>
      <th>log2(# bit ops)</th>
    <!--  <th>log2(# bit ops/it)</th>
      <th>log2(# iterations)</th>-->
      <th>p </th>
      <th>l </th>
      <th>m </th>
      <th>c </th>
      <th>r </th>
      <th>M </th>
      <th># sets A</th>
      </tr>
	  <!-- Berger et al -->
	  <tr>
        <!-- ./isdfq 256 459 255 50 1 3 1 1 1 0 1.00 0 300 -->
        <td>256</td>
        <td>459</td>
        <td>255</td>
        <td>50</td><td> 80 </td>
        <td>77.021000</td>
        <!--<td>21.748776</td>
        <td>55.272224</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <!-- ./isdfq 256 459 255 50 1 3 1 2 1 1 1.30 0 300 -->
        <td>77.099522</td>
        <!--<td>22.726810</td>
        <td>54.372711</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>2</td>
        <td>1</td>
        <td>1.30</td>
        <td>235</td>
      </tr>
	  
      <tr>
        <!-- ./isdfq 256 510 306 50 1 3 1 1 1 0 1.00 0 300 -->
        <td>256</td>
        <td>510</td>
        <td>306</td>
        <td>50</td><td> 90 </td>
        <td>84.786581</td>
        <!--<td>22.034025</td>
        <td>62.752555</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <!-- ./isdfq 256 510 306 50 1 3 1 2 1 1 1.30 0 300 -->
        <td>84.871245</td>
        <!--<td>23.018045</td>
        <td>61.853200</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>2</td>
        <td>1</td>
        <td>1.30</td>
        <td>282</td>
      </tr>
        	  
      <tr>
        <!-- ./isdfq 256 612 408 50 1 3 1 1 1 0 1.00 0 300 -->
        <td>256</td>
        <td>612</td>
        <td>408</td>
        <td>50</td><td> 100 </td>
        <td>98.192084</td>
        <!--<td>22.512521</td>
        <td>75.679563</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <!-- ./isdfq 256 612 408 50 1 3 1 2 1 1 1.30 0 300 -->
        <td>98.288345</td>
        <!--<td>23.507977</td>
        <td>74.780368</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>2</td>
        <td>1</td>
        <td>1.30</td>
        <td>376</td>
      </tr>
        	  
      <tr>
        <!-- ./isdfq 256 765 510 50 1 3 1 1 1 0 1.00 0 300 -->
        <td>256</td>
        <td>765</td>
        <td>510</td>
        <td>50</td><td> 120 </td>
        <td>97.453317</td>
        <!--<td>22.963789</td>
        <td>74.489528</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <!-- ./isdfq 256 765 510 50 1 3 1 2 1 1 1.20 0 300 -->
        <td>97.518929</td>
        <!--<td>23.910554</td>
        <td>73.608374</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>2</td>
        <td>1</td>
        <td>1.20</td>
        <td>433</td>
      </tr>
        	  
      <tr>
        <!-- ./isdfq 1024 450 225 56 1 3 1 1 1 0 1.00 0 300 -->
        <td>1024</td>
        <td>450</td>
        <td>225</td>
        <td>56</td><td> 80 </td>
        <td>76.839227</td>
        <!--<td>23.579323</td>
        <td>53.259904</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <!-- ./isdfq 1024 450 225 56 1 3 1 2 1 1 1.30 0 300 -->
        <td>76.881876</td>
        <!--<td>24.516557</td>
        <td>52.365319</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>2</td>
        <td>1</td>
        <td>1.30</td>
        <td>207</td>
      </tr>
        	  
      <tr>
        <!-- ./isdfq 1024 558 279 63 1 3 1 1 1 0 1.00 0 300 -->
        <td>1024</td>
        <td>558</td>
        <td>279</td>
        <td>63</td><td> 90 </td>
        <td>83.973764</td>
        <!--<td>23.918782</td>
        <td>60.054983</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <!-- ./isdfq 1024 558 279 63 1 3 1 2 1 1 1.20 0 300 -->
        <td>84.001235</td>
        <!--<td>24.798545</td>
        <td>59.202690</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>2</td>
        <td>1</td>
        <td>1.20</td>
        <td>237</td>
      </tr>
        	  
      <tr>
        <!-- ./isdfq 1024 744 372 54 1 3 1 2 1 0 1.00 0 300 -->
        <td>1024</td>
        <td>744</td>
        <td>372</td>
        <td>54</td><td> 110 </td>
        <td>73.645938</td>
        <!--<td>24.940702</td>
        <td>48.705236</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>2</td>
        <td>1</td>
        <!-- ./isdfq 1024 744 372 51 1 3 1 3 1 1 1.30 0 300 -->
        <td>70.542267</td>
        <!--<td>25.640181</td>
        <td>44.902085</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>1.30</td>
        <td>342</td>
      </tr>
        	  
      <!-- Misoczki and Barreto -->
      <tr>
        <!-- ./isdfq 4 2560 1536 128 2 15 1 4 4 0 1.00 0 400 -->
        <td>4</td>
        <td>2560</td>
        <td>1536</td>
        <td>128</td><td> 128 </td>
        <td>181.856969</td>
        <!--<td>27.445102</td>
        <td>154.411867</td>-->
        <td>2</td>
        <td>15</td>
        <td>1</td>
        <td>4</td>
        <td>4</td>
        <!-- ./isdfq 4 2560 1536 128 2 10 1 10 5 1 0.60 0 400 -->
        <td>187.464708</td>
        <!--<td>33.749775</td>
        <td>153.714933</td>-->
        <td>2</td>
        <td>10</td>
        <td>1</td>
        <td>10</td>
        <td>5</td>
        <td>0.60</td>
        <td>288766</td>
      </tr>
        
      <tr>
        <!-- ./isdfq 16 1408 896 128 2 8 1 8 2 0 1.00 0 400 -->
        <td>16</td>
        <td>1408</td>
        <td>896</td>
        <td>128</td><td> 128 </td>
        <td>210.619534</td>
        <!--<td>30.864552</td>
        <td>179.754982</td>-->
        <td>2</td>
        <td>8</td>
        <td>1</td>
        <td>8</td>
        <td>2</td>
        <!-- ./isdfq 16 1408 896 128 2 9 1 10 2 1 1.10 0 400 -->
        <td>210.764955</td>
        <!--<td>31.505106</td>
        <td>179.259849</td>-->
        <td>2</td>
        <td>9</td>
        <td>1</td>
        <td>10</td>
        <td>2</td>
        <td>1.10</td>
        <td>180061</td>
      </tr>
        
      <tr>
        <!-- ./isdfq 256 640 512 64 1 3 1 1 1 0 1.00 0 300 -->
        <td>256</td>
        <td>640</td>
        <td>512</td>
        <td>64</td><td> 102 </td>
        <td>181.673973</td>
        <!--<td>22.849620</td>
        <td>158.824353</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <!-- ./isdfq 256 640 512 64 1 3 1 1 1 1 1.20 0 300 -->
        <td>181.620520</td>
        <!--<td>23.354298</td>
        <td>158.266221</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <td>1.20</td>
        <td>435</td>
      </tr>
        
      <tr>
        <!-- ./isdfq 256 768 512 128 1 3 1 1 1 0 1.00 0 500 -->
        <td>256</td>
        <td>768</td>
        <td>512</td>
        <td>128</td><td> 136 </td>
        <td>253.017954</td>
        <!--<td>23.079168</td>
        <td>229.938785</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <!-- ./isdfq 256 768 512 128 1 3 1 1 1 1 1.20 0 500 -->
        <td>253.010887</td>
        <!--<td>23.633224</td>
        <td>229.377662</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <td>1.20</td>
        <td>435</td>
      </tr>
        
      <tr>
        <!-- ./isdfq 256 1024 512 256 1 3 1 1 1 0 1.00 0 750 -->
        <td>256</td>
        <td>1024</td>
        <td>512</td>
        <td>256</td><td> 168 </td>
        <td>328.997157</td>
        <!--<td>23.451225</td>
        <td>305.545932</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <!-- ./isdfq 256 1024 512 256 1 3 1 1 1 1 1.10 0 750 -->
        <td>329.041324</td>
        <!--<td>23.948907</td>
        <td>305.092416</td>-->
        <td>1</td>
        <td>3</td>
        <td>1</td>
        <td>1</td>
        <td>1</td>
        <td>1.10</td>
        <td>399</td>
      </tr>
      <!-- Misoczki and Barreto -->
      <!-- codes over F2 -->
      <tr>
        <!-- ./isdfq 2 2304 1281 64 2 21 1 8 8 0 1.00 0 300 -->
        <td>2</td>
        <td>2304</td>
        <td>1281</td>
        <td>64</td><td> 80 </td>
        <td>83.556766</td>
        <!--<td>24.383273</td>
        <td>59.173494</td>-->
        <td>2</td>
        <td>21</td>
        <td>1</td>
        <td>8</td>
        <td>8</td>
        <!-- ./isdfq 2 2304 1280 64 2 22 1 9 9 1 0.90 0 300 -->
        <td>83.391304</td>
        <!--<td>24.732941</td>
        <td>58.658363</td>-->
        <td>2</td>
        <td>22</td>
        <td>1</td>
        <td>9</td>
        <td>9</td>
        <td>0.90</td>
        <td>300759</td>
      </tr>

      <tr>
        <!-- ./isdfq 2 3584 1537 128 2 24 1 8 8 0 1.00 0 300 -->
        <td>2</td>
        <td>3584</td>
        <td>1537</td>
        <td>128</td><td> 112 </td>
        <td>112.310859</td>
        <!--<td>24.616879</td>
        <td>87.693980</td>-->
        <td>2</td>
        <td>24</td>
        <td>1</td>
        <td>8</td>
        <td>8</td>
        <!-- ./isdfq 2 3584 1536 128 2 26 1 8 8 1 1.00 0 350 -->
        <td>112.175869</td>
        <!--<td>25.056094</td>
        <td>87.119775</td>-->
        <td>2</td>
        <td>26</td>
        <td>1</td>
        <td>8</td>
        <td>8</td>
        <td>1.00</td>
        <td>481276</td>
      </tr>

      <tr>
        <!-- ./isdfq 2 4096 2048 128 2 25 1 8 8 0 1.00 0 400 -->
        <td>2</td>
        <td>4096</td>
        <td>2048</td>
        <td>128</td><td> 128 </td>
        <td>136.482956</td>
        <!--<td>25.382966</td>
        <td>111.099989</td>-->
        <td>2</td>
        <td>25</td>
        <td>1</td>
        <td>8</td>
        <td>8</td>
        <!-- ./isdfq 2 4096 2048 128 2 26 1 9 9 1 1.00 0 400 -->
        <td>136.470143</td>
        <!--<td>25.973548</td>
        <td>110.496595</td>-->
        <td>2</td>
        <td>26</td>
        <td>1</td>
        <td>9</td>
        <td>9</td>
        <td>1.00</td>
        <td>855741</td>
      </tr>

      <tr>
        <!-- ./isdfq 2 7168 3073 256 3 37 1 32 8 0 1.00 0 600 -->
        <td>2</td>
        <td>7168</td>
        <td>3073</td>
        <td>256</td><td> 192 </td>
        <td>216.068572</td>
        <!--<td>35.617758</td>
        <td>180.450814</td>-->
        <td>3</td>
        <td>37</td>
        <td>1</td>
        <td>32</td>
        <td>8</td>
        <!-- ./isdfq 2 7168 3072 256 3 38 1 36 9 1 1.00 0 700 -->
        <td>215.907793</td>
        <!--<td>36.465708</td>
        <td>179.442085</td>-->
        <td>3</td>
        <td>38</td>
        <td>1</td>
        <td>36</td>
        <td>9</td>
        <td>1.00</td>
        <td>1079376989</td>
      </tr>

      <tr>
        <!-- ./isdfq 2 8192 4097 256 3 38 1 16 8 0 1.00 0 700 -->
        <td>2</td>
        <td>8192</td>
        <td>4097</td>
        <td>256</td><td> 256 </td>
        <td>265.158785</td>
        <!--<td>36.933233</td>
        <td>228.225551</td>-->
        <td>3</td>
        <td>38</td>
        <td>1</td>
        <td>16</td>
        <td>8</td>
        <!-- ./isdfq 2 8192 4096 256 3 39 1 17 9 1 1.00 0 700 -->
        <td>265.007138</td>
        <!--<td>37.777370</td>
        <td>227.229767</td>-->
        <td>3</td>
        <td>39</td>
        <td>1</td>
        <td>17</td>
        <td>9</td>
        <td>1.00</td>
 	<!-- here (k,p)/sqrt((2p,p)) is bigger than 2^32 -->
	<td>2561023536</td>
      </tr>
    </table>
    (the commands to verify these parameters are given as comments in the source code of this html document)
  </li>
</ul>

<p><a name="m">(*)</a>
The algorithm looks for words having weight <tt>p</tt> among the
positions indexed by a set X, weight <tt>p</tt> among the positions indexed by
a set Y and weight 0 among <tt>l</tt> positions indexed by a set Z. <br></br>

In the paper <a href="./#2008.mceliece">Attacking and defending the McEliece
cryptosystem</a> Dan Bernstein, Tanja Lange, and I 
proposed to speed up Stern's algorithm for codes over F<sub>2</sub> by
taking <tt>m</tt> sets Z<sub>1</sub>,...,Z<sub>m</sub> outside the
information set, each with <tt>l</tt> positions
and to look for words having weight <tt>p</tt> among the X-indexed positions, weight <tt>p</tt> among the Y-indexed positions, and weight 0 among the positions indexed by one
of those sets Z<sub>j</sub>. For
small q, such as q=2, it turns out that using <tt>m>1</tt> provides better
results.
<br></br>

For codes over arbitrary fields F<sub>q</sub> I considered only one set Z. That's why
the default setting of 
<a href="./isdfq.c">isdfq.c</a> is <tt>m=1</tt>.
</p>


<!-- Last updated -->
<p class="right">
<br></br>
<script type="text/javascript">
var LastUpdated = document.lastModified;
document.writeln ("This website was last updated on " + LastUpdated);
</script>
</p> 

</body>

</html>

