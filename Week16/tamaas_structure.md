Since we are contributing [Tammas](https://gitlab.com/tamaas/tamaas), we first need to figure out the structure of tamaas works.

### src

#### model

##### model class

Model containing pressure and displacement. This class is a container for the model fields. It is supposed to be dimension agnostic, hence the GridBase members.

`model.hh` contains getDiscretization() function, we can use this function for *n* and *m* variable of our python script.


