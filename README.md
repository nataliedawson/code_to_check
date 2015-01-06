OrengoGroup Code Review 1
=========================

This repo provides access to a Perl script that has been submitted for review within the Orengo group. Instructions for testing on your local system are given below. 

Getting a local copy of the code
--------------------------------

Make sure you've got git installed (you probably don't need to do this)

    $ sudo apt-get install git

Create a sensible place in your home directory

    $ mkdir ~/orengogroup_code_review
    $ cd ~/orengogroup_code_review
  
Clone this repository (i.e. get a local copy that you can play with):

    $ git clone https://github.com/nataliedawson/code_to_check.git
    Cloning into 'code_to_check'...
    remote: Counting objects: 10, done.
    remote: Compressing objects: 100% (8/8), done.
    remote: Total 10 (delta 0), reused 6 (delta 0)
    Unpacking objects: 100% (10/10), done.
    Checking connectivity... done.

Move into this directory:

    $ cd code_to_check

Using the right Perl
--------------------

This Perl script relies on a bunch of libraries being available (both from CPAN and local CATH libraries). Luckily we maintain a version of Perl in a shared location that already has these requirements sorted out for you.

Find out which perl you are running.

    $ which perl
    /opt/local/perls/build-trunk/bin/perl

If you get a different answer to above then you can alter your `PATH` to make the shared perl available by default.

    # in bash
    $ PATH=/opt/local/perls/build-trunk/bin:${PATH}

If you have problems with any of the above then speak to someone who uses perl in the group (e.g. Natalie, Ian, Tony, etc).

Running the script
------------------

    $ perl get_rep_for_all_families.pl -i data/115699.faa -o data/ -d data/ -v 4.0

What the code does
------------------

The algorithm aims to identify the most appropriate representative for a given FunFam alignment.

If the FunFam contains one or more CATH domains:

  * identify the CATH domain with the highest accumulated structural similarity (SSAP) score

Otherwise:

  * identify the Gene3D sequence with the highest (???) GO score (not yet fully implemented)

What I would like people to check
---------------------------------

 1. issue 1
 1. issue 2
 1. issue 3


