==================
HELICITY ANALYZER
==================

---------------------
Author: C. Yero
Date: Oct 18, 2020

email: cyero@jlab.org | cyero002@fiu.edu
---------------------

Brief: The Helicity Analyzer, helicityAnalyzer.cpp(.h),
is a derived C++ class from the more generic baseAnalyzer.cpp(.h)
The idea is that the baseAnalyzer serves as a general purpose analyzer
which does essential tasks that every experimental analysis in Hall C
should do.

The more specific (derived) classes such as the helicity class,
have specific aanalysis requirement, but also can carry out the
general analysis requiremts. One can work on specific methods
of a derived class without altering the base class.

If one changes the base class, however, it affects all classes that
derive from it. This generalization of having a single base class is
preferred that having to copy/paste specific functions, which is prone
to making mistakes.

----------------------

I will keep updating the README files / Documentation as
progress is made.