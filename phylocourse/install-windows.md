# Install and run FastME on Windows

(1) **Obtaining FastME.** Go to [http://www.atgc-montpellier.fr/fastme/binaries.php](http://www.atgc-montpellier.fr/fastme/binaries.php).
At the bottom of the page you'll find a Download button (you don't have to fill
in your name and email). Hit download and save the associated
`fastme-2.1.5.tar.gz` file to some directory (presumable your `Downloads`
folder).


(2) **Install FastME.** Make a new folder somewhere in your filesystem where you
find it most convenient. I made a folder `fastme-exercise`. Copy the `.tar.gz` file in the new directory. Now unzip the `.tar.gz` archive. This can be done in two ways:

1. Either use a free software like [`7-zip`](https://www.7-zip.org/) or
2. Use the command prompt. Since we will be operating from the command line
   later, it might be worthwhile to try this out. On relatively recent versions
   of Windows (I guess since Win 10?), you can apparently use Linux style commands in the command prompt. To do so, open the file explorer in the folder where you put the `.tar.gz` file. Then open the command prompt by holding the shift key and right-clicking in the file manager window, you'll see something like this:

   ![](/assets/phylocourse/bbwin/1.png)

   Click on the `Open powershell window here` option. A scary black screen will
   appear like so:

   ![](/assets/phylocourse/bbwin/2.png)

   Now type `tar -xvf fastme-2.1.5.tar.gz` (tip: typing `tar -xvf fastme` and
   hitting the tab key will provide some hints [so-called 'tab-completion']).
   Hit `Enter` and the `tar.gz` archive will unpack.

   ![](/assets/phylocourse/bbwin/3.png)

After you successfully unpacked the `.tar.gz` archive, you can copy the executable
`fastme.exe` in the `binaries` folder of the `fastme-2.1.5` folder in the directory
you made. You should now have a view like this in your file manager:

![](/assets/phylocourse/bbwin/4.png)

(3) **Run FastME.** To run FastME on some data, put the data in the same directory
as where the `fastme.exe` executable is located. For instance, you can copy the
`nucleic.phy` file from the `examples` folder in the `fastme-2.5.1` directory
to your directory where the `fastme.exe` file is located:

![](/assets/phylocourse/bbwin/5.png)

Now open the Powershell in the same folder (so the one where the `fastme.exe`
file and data file(s) are located) again like explained above (`Shift + right
click` in the file manager). Now type `fastme.exe -i nucleic.phy -dJC69 -m NJ` to
run tree inference using Neighbor-Joining with Jukes & Cantor based distances

![](/assets/phylocourse/bbwin/6.png)

Now you have generated some output files, which you can examine in the file
manager:

![](/assets/phylocourse/bbwin/7.png)

# Installing and running IQ-TREE on Windows

I would recommend the exact same approach as above, but now download IQ-TREE
[here](http://www.iqtree.org/#download). After unzipping, you'll find the
`iqtree.exe` file in the `bin` folder of the `iqtree-<version>-Windows`
folder.
