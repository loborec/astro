php-astro
===========

Class Astro is small collection of the astronomical functions written in php. 

This class is based on the works:

Darin Koblick https://www.mathworks.com/matlabcentral/profile/authors/1284781-darin-koblick

Dubravko Loborec http://www.dubravkoloborec.com


## Installing

You can use [Composer](http://getcomposer.org/) to add the [package](https://packagist.org/packages/loborec/zf2-message) to your project:

```json
{
  "require": {
    "loborec/astro": "dev-master"
  }
}
```

## Functions
# Lunar Azimuth and Altitude Estimation Algorithm (LunarAzEl)
Predict the Lunar Azimuth and Altitude within +/- .2 deg of any lat and lon for a given UTC

This algorithm will accept a Latitude, Longitude and Altitude location as well as a specific universal coordinated time. It will use this information and calculate the position of the moon in a local coordinate frame (az and alt aka az and el).

Example Function Call: 
```php
$result = Astro::LunarAzEl(1991, 05, 19, 13, 0, 0, 50, 10, 0);
var_dump(result);
```
Result:
```txt
array (size=2)
  0 => float 110.63005216655
  1 => float 35.226177186756
```

Input List: 

UTC Date and Time - Year, Month, Day, Hour, Minute, Second

Latitude - Site Latitude in degrees -90:90 -> S(-) N(+) 

Longitude - Site Longitude in degrees -180:180 W(-) E(+) 

Altitude - Site Altitude in km

Output List: 

Az - Lunar Azimuth angle in degrees 

El - Lunar Elevation/Altitude Angle in degrees

Verified output by comparison with the following source data: http://aa.usno.navy.mil/data/docs/AltAz.php

Copyright (c) 2010, Darin Koblick


## Examples



            



