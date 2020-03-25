EXTERNAL SOFTWARE CONNECTION
============================

A folder for each connector.

The importance of the **configSBI.txt** file:
If the parameters of the program instalation are others, a configSBI file must be created.

Then it can be added to your scripts by:

<code>
import os

os.environ['SBI_CONFIG_FILE'] = "your_file"
</code>