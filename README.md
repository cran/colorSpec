## colorSpec


**colorSpec** is an **R** package providing an S3 class with methods for color spectra. 
It supports the standard calculations with spectral properties of light sources, materials, cameras, eyes, scanners, etc. 
And it works well with the more general action spectra.



When a spectrum is contructed, a `quantity` is required,
though it can be changed after contruction.
The **colorSpec** quantities are typically not the same as the SI quantities;
they are more general.
There are 14 quantities, which then defined the 4 basic spectrum `types`,
as given in this table from the User Guide.

<TABLE WIDTH=856 BORDER=2 BORDERCOLOR="#000000" CELLPADDING=5 CELLSPACING=0>
<COL WIDTH=116>
<COL WIDTH=163>
<COL WIDTH=91>
<COL WIDTH=227>
<COL WIDTH=205>
<TR VALIGN=BOTTOM>
<TH WIDTH=116>
<P><FONT FACE="Arial, sans-serif"><FONT SIZE=2>colorSpec<BR><FONT FACE="Courier New, monospace">type</FONT></FONT></FONT></P>
</TH>
<TH WIDTH=163>
<P><FONT FACE="Courier New, monospace"><FONT SIZE=2>quantity</FONT></FONT></P>
</TH>
<TH WIDTH=91>
<P><FONT FACE="Arial, sans-serif"><FONT SIZE=2>metric</FONT></FONT></P>
</TH>
<TH WIDTH=227>
<P><FONT FACE="Arial, sans-serif"><FONT SIZE=2>comments</FONT></FONT></P>
</TH>
<TH WIDTH=205>
<P><FONT FACE="Arial, sans-serif"><FONT SIZE=2>examples<BR>(objects,
files, functions)</FONT></FONT></P>
</TH>
</TR>
<TR>
<TD ROWSPAN=2 WIDTH=116>
<P ALIGN=CENTER><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">light</FONT></FONT></P>
</TD>
<TD WIDTH=163>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">energy</FONT></FONT></P>
</TD>
<TD WIDTH=91>
<P ALIGN=CENTER STYLE="margin-top: 0.1in"><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">radiometric</FONT></FONT></P>
</TD>
<TD WIDTH=227>
<P ALIGN=LEFT><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">radiometric
quantities are conventional in colorimetry</FONT></FONT></P>
</TD>
<TD WIDTH=205>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">D65.1nm<BR>pos1-20x.scope<BR>BlueFlame.txt</FONT></FONT></P>
</TD>
</TR>
<TR>
<TD WIDTH=163>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">photons</FONT></FONT></P>
</TD>
<TD WIDTH=91>
<P ALIGN=CENTER STYLE="margin-top: 0.1in"><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">actinometric</FONT></FONT></P>
</TD>
<TD WIDTH=227>
<P ALIGN=LEFT><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">for
color calculations, actinometric units are automatically
converted to radiometric</FONT></FONT></P>
</TD>
<TD WIDTH=205>
<FONT FACE="Courier New"><FONT SIZE=1 STYLE="font-size: 8pt">F96T12</FONT></FONT><BR>
<FONT FACE="Courier New"><FONT SIZE=1 STYLE="font-size: 8pt">Airam-GR8E.txt</FONT></FONT>
</TD>
</TR>
<TR>
<TD ROWSPAN=6 WIDTH=116>
<P ALIGN=CENTER><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">responsivity<BR>.light</FONT></FONT></P>
</TD>
<TD WIDTH=163>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">energy-&gt;electrical</FONT></FONT></P>
</TD>
<TD ROWSPAN=3 WIDTH=91>
<P ALIGN=CENTER STYLE="margin-top: 0.1in"><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">radiometric</FONT></FONT></P>
</TD>
<TD WIDTH=227>
<P ALIGN=LEFT><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">RGB
camera response</FONT></FONT></P>
</TD>
<TD WIDTH=205>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">Flea2.RGB<BR>Red-Epic-Dragon.txt</FONT></FONT></P>
</TD>
</TR>
<TR>
<TD WIDTH=163>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">energy-&gt;neural</FONT></FONT></P>
</TD>
<TD WIDTH=227>
<P ALIGN=LEFT><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">eye
response</FONT></FONT></P>
</TD>
<TD WIDTH=205>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">xyz1931.1nm<BR>luminsivity.1nm<BR>Osmia-rufa.txt</FONT></FONT></P>
</TD>
</TR>
<TR>
<TD WIDTH=163>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">energy-&gt;action</FONT></FONT></P>
</TD>
<TD WIDTH=227>
<P ALIGN=LEFT><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">examples
are erythemal action, melatonin suppression, etc.</FONT></FONT></P>
</TD>
<TD WIDTH=205>
<P ALIGN=LEFT STYLE="font-style: normal; font-weight: normal"><FONT COLOR="#000000"><FONT FACE="Courier New"><FONT SIZE=1 STYLE="font-size: 8pt">erythemalSpectrum()</FONT></FONT></FONT></P>
</TD>
</TR>
<TR>
<TD WIDTH=163>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">photons-&gt;electrical</FONT></FONT></P>
</TD>
<TD ROWSPAN=3 WIDTH=91>
<P ALIGN=CENTER STYLE="margin-top: 0.1in"><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">actinometric</FONT></FONT></P>
</TD>
<TD WIDTH=227>
<P ALIGN=LEFT><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">silicon
sensors usually use quantum efficiency as measure of responsivity</FONT></FONT></P>
</TD>
<TD WIDTH=205>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">Zyla_sCMOS.txt<BR>FoveonX3.txt</FONT></FONT></P>
</TD>
</TR>
<TR>
<TD WIDTH=163>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">photons-&gt;neural</FONT></FONT></P>
</TD>
<TD WIDTH=227>
<P ALIGN=LEFT><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">response
units might be photocurrent, or spikes/sec, etc.</FONT></FONT></P>
</TD>
<TD WIDTH=205>
<FONT FACE="Courier New"><FONT SIZE=1 STYLE="font-size: 8pt">HigherPasserines</FONT></FONT>
</TD>
</TR>
<TR>
<TD WIDTH=163>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">photons-&gt;action</FONT></FONT></P>
</TD>
<TD WIDTH=227>
<P ALIGN=LEFT><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">photosynthesis
is an example</FONT></FONT></P>
</TD>
<TD WIDTH=205>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">BeanPhotosynthesis.txt</FONT></FONT></P>
</TD>
</TR>
<TR>
<TD ROWSPAN=2 WIDTH=116>
<P ALIGN=CENTER><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">material</FONT></FONT></P>
</TD>
<TD WIDTH=163>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">reflectance</FONT></FONT></P>
</TD>
<TD WIDTH=91>
<P ALIGN=CENTER STYLE="margin-top: 0.1in"><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">NA</FONT></FONT></P>
</TD>
<TD WIDTH=227>
<P ALIGN=LEFT><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">this
is the diffuse reflectance</FONT></FONT></P>
</TD>
<TD WIDTH=205>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">CC_Avg20_spectrum_XYY.txt</FONT></FONT></P>
</TD>
</TR>
<TR>
<TD WIDTH=163>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">transmittance<BR>absorbance</FONT></FONT></P>
</TD>
<TD WIDTH=91>
<P ALIGN=CENTER STYLE="margin-top: 0.1in"><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">NA</FONT></FONT></P>
</TD>
<TD WIDTH=227>
<P ALIGN=LEFT><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">for
color calculations, <I>absorbance</I> is automatically converted
to <I>transmittance</I></FONT></FONT></P>
</TD>
<TD WIDTH=205>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">Hoya<BR>Hematoxylin.txt<BR>atmosphere2003</FONT></FONT></P>
</TD>
</TR>
<TR>
<TD WIDTH=116>
<P ALIGN=CENTER><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">responsivity<BR>.material</FONT></FONT></P>
</TD>
<TD WIDTH=163>
<P ALIGN=LEFT><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">material-&gt;electrical<BR>material-&gt;neural<BR>material-&gt;action</FONT></FONT></P>
</TD>
<TD WIDTH=91>
<P ALIGN=CENTER STYLE="margin-top: 0.1in"><FONT FACE="Courier New, monospace"><FONT SIZE=1 STYLE="font-size: 8pt">NA</FONT></FONT></P>
</TD>
<TD WIDTH=227>
<P ALIGN=LEFT><FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">a
spectrum of this type typically comes from both a light source
and a camera</FONT></FONT></P>
</TD>
<TD WIDTH=205>
<FONT COLOR="#000000"><FONT FACE="Courier New"><FONT SIZE=1 STYLE="font-size: 8pt">scanner.ACES</FONT></FONT></FONT>
<br>
<FONT FACE="Arial, sans-serif"><FONT SIZE=1 STYLE="font-size: 8pt">(a standard for scanning film)</FONT></FONT>
</TD>
</TR>
</TABLE>

<P ALIGN=CENTER STYLE="margin-top: 0.08in; margin-bottom: 0in"><FONT FACE="Arial, sans-serif"><FONT SIZE=2><FONT SIZE=3>
Table 2.1.  The </FONT><FONT FACE="Courier New, monospace"><FONT SIZE=3>types</FONT></FONT>
<FONT SIZE=3>of spectra and their </FONT><FONT FACE="Courier New, monospace"><FONT SIZE=3>quantities</FONT></FONT>
</FONT></FONT>
</P>

<br><br>
The `quantity` is used to label plots,
and to make sense out of the arguments in the function `product()`.
These quantities are typically not the same as the SI quantities; they are more general.
For details please see the **colorSpec User Guide**.

<br><br>
The function `invert()` implements a method in 
<blockquote>
<p>Davis G (2019).
&ldquo;A Centroid for Sections of a Cube in a Function Space, with application to Colorimetry.&rdquo;
<em>ArXiv e-prints</em>.
1811.00990, <a href="https://arxiv.org/abs/1811.00990">https://arxiv.org/abs/1811.00990</a>. 
</p>
</blockquote>

