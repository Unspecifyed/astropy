from astropy.coordinates import get_body, EarthLocation, AltAz, get_sun, SkyCoord
from astropy.time import Time
import datetime
import astropy.units as u
from astropy.coordinates import solar_system_ephemeris, get_body
from astropy.coordinates import ICRS, FK5, Galactic, GeocentricTrueEcliptic
from astropy.coordinates import Angle

# Function to determine which zodiac constellation a given position is in, based on ICRS coordinates
def get_zodiac_sign_icrs(ecliptic_longitude):
    zodiac_boundaries = [
        ('Pisces', 29.0),
        ('Aries', 53.0),
        ('Taurus', 90.0),
        ('Gemini', 118.0),
        ('Cancer', 138.0),
        ('Leo', 173.0),
        ('Virgo', 217.0),
        ('Libra', 242.0),
        ('Scorpio', 247.0),
        ('Ophiuchus', 266.0),
        ('Sagittarius', 295.0),
        ('Capricorn', 326.0),
        ('Aquarius', 352.0),
        ('Pisces', 360.0)
    ]

    for i in range(len(zodiac_boundaries) - 1):
        if zodiac_boundaries[i][1] <= ecliptic_longitude < zodiac_boundaries[i + 1][1]:
            return zodiac_boundaries[i][0]
    if ecliptic_longitude >= zodiac_boundaries[-1][1] or ecliptic_longitude < zodiac_boundaries[0][1]:
        return zodiac_boundaries[0][0]
    return 'Unknown'

# Function to determine the current season and time until the next equinox/solstice
def get_season_info(now):
    year = now.datetime.year
    seasons = {
        'spring_equinox': Time(f'{year}-03-20T00:00:00', format='isot', scale='utc'),
        'summer_solstice': Time(f'{year}-06-21T00:00:00', format='isot', scale='utc'),
        'fall_equinox': Time(f'{year}-09-23T00:00:00', format='isot', scale='utc'),
        'winter_solstice': Time(f'{year}-12-21T00:00:00', format='isot', scale='utc')
    }

    if seasons['spring_equinox'] <= now < seasons['summer_solstice']:
        current_season = 'Spring'
        next_event = 'Summer Solstice'
        time_until_next = seasons['summer_solstice'] - now
    elif seasons['summer_solstice'] <= now < seasons['fall_equinox']:
        current_season = 'Summer'
        next_event = 'Fall Equinox'
        time_until_next = seasons['fall_equinox'] - now
    elif seasons['fall_equinox'] <= now < seasons['winter_solstice']:
        current_season = 'Fall'
        next_event = 'Winter Solstice'
        time_until_next = seasons['winter_solstice'] - now
    else:
        current_season = 'Winter'
        next_event = 'Spring Equinox'
        time_until_next = seasons['spring_equinox'].replace(year=year + 1) - now

    return current_season, next_event, time_until_next

# Set observer location
def track_planets():
    location = EarthLocation(lat=51.5074*u.deg, lon=-0.1278*u.deg, height=0*u.m)  # Change lat/lon for your location
    now = Time(datetime.datetime.now(datetime.UTC))
    frame = AltAz(obstime=now, location=location)
    planets = ['mercury', 'venus', 'mars', 'jupiter', 'saturn', 'uranus']  # 'pluto' is commented out
    additional_bodies = ['sun', 'moon']

    current_season, next_event, time_until_next = get_season_info(now)
    print(f"Current season: {current_season}")
    print(f"Next event: {next_event} in {time_until_next}")

    with solar_system_ephemeris.set('builtin'):
        for planet in planets:
            try:
                body = get_body(planet, now, location)
            except KeyError as e:
                if 'pluto' in str(e):
                    print("Pluto's position cannot be calculated with the 'builtin' ephemeris.")
                    continue
            body_altaz = body.transform_to(frame)
            ecliptic = body.transform_to(GeocentricTrueEcliptic())
            ecliptic_longitude = ecliptic.lon.deg
            zodiac_sign = get_zodiac_sign_icrs(ecliptic_longitude)

            print(f"{planet.capitalize()} is {body_altaz.alt:.2f} degrees above the horizon, located in the zodiac constellation: {zodiac_sign}")

        for body_name in additional_bodies:
            if body_name == 'sun':
                body = get_sun(now)
            elif body_name == 'moon':
                body = get_body('moon', now, location)
            
            body_altaz = body.transform_to(frame)
            ecliptic = body.transform_to(GeocentricTrueEcliptic())
            ecliptic_longitude = ecliptic.lon.deg
            zodiac_sign = get_zodiac_sign_icrs(ecliptic_longitude)

            print(f"{body_name.capitalize()} is {body_altaz.alt:.2f} degrees above the horizon, located in the zodiac constellation: {zodiac_sign}")

if __name__ == "__main__":
    track_planets()
