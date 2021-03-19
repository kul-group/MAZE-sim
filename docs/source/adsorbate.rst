Adsorbate Class
================

The Adsorbate module introduces a class specific to adsorbates which can be added and removed from the host zeolite.

Creating an Adsorbate Object
------------------------------

For this example, we will begin with a BEA zeolite with a Sn atom.

 .. code:: python

    # Create a zeotype from the BEA structure
    iz = ImpefectZeotype.make('BEA')
    # Change atom 186 to a Sn atom
    iz[186].symbol = 'Sn'

Now we can create an adsorbate object, in this case using an ammonia molecule:

 .. code:: python

    mol = molecule('NH3')
    ads = Adsorbate(mol, host_zeotype=iz)

Positioning Adsorbates
----------------------------------

So far, we have created a Zeotype object with a BEA structure as well as an Adsorbate object with an ammonia structure. Now, the ammonia molecule can be added to the BEA framework. The `position_ads()` function can take user-specified host binding site, adsorbate binding location, and adsorbate binding atom to position and orient the adsorbate.

 .. code:: python

    ads2 = ads.position_ads(donor_ind=0, host_ind=186, pos=[5.6, 7.8, 15.0]

Specifying all of these parameters is often time-consuming. The `position_ads()` function has the ability to automatically select all of these parameters.

 .. code:: python

    ads2 = ads.position_ads()

It is important to note that the automated placement method merely provides a good guess for adsorbate placement. For larger adsorbates, this placement may be incorrect, so these structures should be checked before use.


Adding and Removing Adsorbates
-----------------------------------

With our adsorbate molecule properly positioned, we can now add it to the zeolite framework using the `integrate_adsorbate()` function.

 .. code:: python

    iz.integrate_adsorbate(ads2)

The adsorbate can be removed from the zeolite framework using the `remove_adsorbate()` function:


 .. code:: python

    zeotype.remove_adsorbate(ads2)

Note that multiple adsorbates can be added and removed independently of one another.