.. _potable-troubleshooting:

Troubleshooting Potable Input Files
===================================

Checking tabulated functions
----------------------------

Sometimes you will want to check if the tabulated functions defining your potential model are as they should be. Perhaps the easiest way to do this is to use the ``excel`` tabulation targets:

    * **excel:** dumps pair potential models into an Excel spreadsheet file,
    * **excel_eam:** as above but for EAM potentials,
    * **excel_eam_fs:** for Finnis-Sinclair EAM potential models.

If given as the value of the ``target`` field in a potable input file's ``[Tabulation]`` section each one will dump the model into an Excel ``.xlsx`` formatted spreadsheet. 

Depending on the type of model the spreadsheet can contain **Pair**\ , **EAM-Density** and **EAM-Embed** sheets. The first column of each is either the *r* or *rho* values of the function with the remaining columns giving the function values.

To enable quick checks it is often more convenient to use potable's ``--override-item`` option to temporarily set the target, rather than having to edit the input file. For example, it is used here to set the ``excel`` target::

    potable --override-item=Tabulation:target=excel potential_model.aspot test.xlsx



