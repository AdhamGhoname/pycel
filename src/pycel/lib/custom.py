import math
from pycel.excellib import _numerics, DIV0, NUM_ERROR
import numpy as np

c1 = 6.112
c2 = 17.67
c3 = 243.5
c4 = 0.00066
c5 = 0.000000759
c8 = -10440.4
c9 = -11.29465
c10 = -0.027022355
c11 = 0.00001289036
c12 = -0.000000002478068
c13 = 6.5459673
c14 = 100.45
c15 = 33.193
c16 = 2.319
c17 = 0.17074
c18 = 1.2063
e = 2.71828


def get_wsat(wb, elevInFt):
    try:
        wb = float(wb)
        elevInFt = float(elevInFt)
        # Wet bulb in Degree Rankin
        rt = wb + 459.67

        # Atmos. Pressure in psia
        pt = 14.696 * (1 - 0.0000068753 * elevInFt) ** 5.2559

        pws = math.exp(
            c8 / rt + c9 + c10 * rt + c11 * rt ** 2 + c12 * rt ** 3 + c13 * math.log(rt)
        )

        wsat = (pws * 0.62198) / (pt - pws)
    except :
        return DIV0, DIV0
    return pt, wsat


def dewpoint(db, rh):
    try:
        db = float(db)
        rh = float(rh)
        es = 6.11 * 10 ** ((7.5 * ((db - 32) * 5 / 9)) / (237.7 + ((db - 32) * 5 / 9)))
        x = math.log(es * rh / 611)
        yd = math.log(10)
        zd = x / yd
        dw = (237.7 * zd) / (7.5 - zd)
        dw = (dw * 9 / 5) + 32
    except :
        return DIV0
    return dw


def dewpointal(elev, hr):
    try:
        elev = float(elev)
        hr = float(hr)
        p = 14.696 * (1 - 0.0000068753 * elev) ** 5.2559
        pw = (p * hr) / (0.62198 + hr)
        alpha = math.log(pw)
        final = c14 + c15 * alpha + c16 * alpha ** 2 + c17 * alpha ** 3 + c18 * (pw) ** 0.1984
    except :
        return DIV0
    return final


def relhum(db, wb, elevInFt):
    try:
        db = float(db)
        wb = float(wb)
        elevInFt = float(elevInFt)
        pt, wsat = get_wsat(wb, elevInFt)
        if DIV0 in [pt, wsat]:
            return NUM_ERROR
        wnom = (1093 - 0.556 * wb) * wsat - 0.24 * (db - wb)
        wdenom = 1093 + 0.444 * db - wb
        w = wnom / wdenom
        pw = (w * pt) / (0.62198 + w)
        rt = db + 459.67
        ps = math.exp(
            c8 / rt
            + c9
            + c8 * rt
            + c10 * rt
            ** 2 + c11 * rt
            ** 3 + c12 * math.log(rt)
        )
        pa = pt - pw
        return 100 * (w * pa) / (0.62198 * ps)
    except :
        return DIV0


def humrat(db, wb, elevInFt):
    try:
        db = float(db)
        wb = float(wb)
        elevInFt = float(elevInFt)

        pt, wsat = get_wsat(wb, elevInFt)
        if DIV0 in [pt, wsat]:
            return NUM_ERROR
        wnom = (1093 - 0.556 * wb) * wsat - 0.24 * (db - wb)
        wdenom = 1093 + 0.444 * db - wb
        return wnom / wdenom
    except :
        return DIV0


def dryfromhr(hr, wb, elevInFt):
    try:
        hr = float(hr)
        wb = float(wb)
        elevInFt = float(elevInFt)
        pt, wsat = get_wsat(wb, elevInFt)
        if DIV0 in [pt, wsat]:
            return NUM_ERROR
        a1 = (1093 - 0.556 * wb) * wsat
        a2 = 1093 - wb
        a3 = hr * 0.444
        db = (a1 - hr * a2 + 0.24 * wb) / (0.24 + a3)
    except :
        return DIV0
    return db


def maxwetc(dbr, dbmax, eff):
    try:
        dbr = float(dbr)
        dbmax = float(dbmax)
        eff = float(eff)
        wmax = dbmax
        for _ in range(1000):
            dbi = dbmax - ((dbmax - wmax) * eff)
            wmax = wmax - 0.1
            if dbi < dbr:
                break
    except :
        return DIV0
    return wmax


def loadcalc(db, dp, cfm, win):
    try:
        db = float(db)
        dp = float(dp)
        cfm = float(cfm)
        win = float(win)
        """Calc for the coil 6 Row 500FPM load"""
        fc = 0.38
        gpm = (cfm / 2500) * 12.5
        wout1 = (db - win) * fc
        wout2 = fc * (dp - win) if win < dp else 0
        wout = win + wout1 + wout2
        btu = (wout - win) * gpm * 8 * 60
        ton = btu / 12000
    except :
        return DIV0
    return ton


def st(hr, elev, stemp, eh):
    try:
        hr = float(hr)
        elev = float(elev)
        stemp = float(stemp)
        eh = float(eh)

        hr1 = hr
        stemp1 = stemp
        e1 = 0.24 * stemp + hr * (1061 + 0.444 * stemp)
        while e1 <= eh:
            stemp = stemp - 0.001
            hr = humrat(stemp, stemp, elev)
            if hr == DIV0:
                return NUM_ERROR
            e1 = 0.24 * stemp + hr * (1061 + 0.444 * stemp)
        st = round(stemp, 3)
        hr = None
    except :
        return DIV0
    return stemp1


def heatindex(db, rh):
    try:
        db = float(db)
        rh = float(rh)

        if db <= 40:
            heatindex = db

        if db > 40:
            heatindex = (-10.3) + (1.1 * db) + (0.047 * rh)

        if db >= 79:
            heatindex = (
                -42.379
                + (2.04901523 * db)
                + (10.14333127 * rh)
                - (0.22475541 * db * rh)
                - (6.83783 * (10 ** -3) * (db ** 2))
                - (5.481717 * (10 ** -2) * (rh ** 2))
                + (1.22874 * (10 ** -3) * (db ** 2) * rh)
                + (8.5282 * (10 ** -4) * db * (rh ** 2))
                - (1.99 * (10 ** -6) * (db ** 2) * (rh ** 2))
            )

        j1 = ((13 - rh) / 4) * ((17 - abs(db - 95)) / 17)
        j2 = ((rh - 85) / 10) * ((87 - db) / 5)

        cn1 = 0
        cn2 = 0

        if db >= 80:
            cn1 = cn1 + 1
            cn2 = cn2 + 1
        if db <= 112:
            cn1 = cn1 + 1
        if rh < 13:
            cn1 = cn1 + 1
        if cn1 >= 3:
            heatindex = heatindex - j1

        if db <= 87:
            cn2 = cn2 + 1
        if rh > 85:
            cn2 = cn2 + 1
        if cn2 >= 3:
            heatindex = heatindex + j2
    except :
        return DIV0
    return heatindex


def endtemp(db, wb, elev, roomload, cfm, dxextra, deff, Ideff, suptemp):
    try:
        db = float(db)
        wb = float(wb)
        elev = float(elev)
        roomload = float(roomload)
        cfm = float(cfm)
        dxextra = float(dxextra)
        deff = float(deff)
        Ideff = float(Ideff)
        suptemp = float(suptemp)

        dbo = db

        if db > suptemp:
            if (db - (db - wb) * deff) > suptemp:
                if (db - (db - wb) * Ideff) > suptemp:
                    senfactoer = 1.08 - (3.45454545454546e-05 * elev)
                    latentfactoer = 4.45 - (0.00014234 * elev)
                    # start up calculate how much idec can supply
                    rhe = relhum(db, wb, elev)
                    Dp = dewpoint(db, rhe)
                    hre = humrat(db, wb, elev)
                    eh = enthlpywet(db, wb, elev)
                    fstload = loadcalc(db, Dp, cfm, wb + 1) + dxextra
                    ehr = eh - ((fstload * 12000) / (latentfactoer * cfm))
                    maxsen = (senfactoer * cfm * (db - Dp)) / 12000
                    sttemp = st(hre, elev, Dp, ehr)
                    senstemp = db - ((fstload * 12000) / (senfactoer * cfm))
                    if maxsen >= fstload:
                        finltemp = senstemp
                    else:
                        finltemp = sttemp
                    finlwet = weten(ehr, elev)
                    finleh = ehr
                    finldry = finltemp - (finltemp - finlwet) * deff
                    roomeh = finleh + ((roomload * 12000) / (latentfactoer * cfm))
                    roomwet = weten(roomeh, elev)
                    roomdry = ((roomload * 0.7 * 12000) / (senfactoer * cfm)) + finldry
                    roomeh1 = roomeh
                    roomwet1 = roomwet
                    roomdry1 = roomdry
                    if finldry > suptemp:
                        for _ in range(1000):
                            roomeh2 = roomeh1
                            roomwet2 = roomwet1
                            roomdry2 = roomdry1
                            secload = loadcalc(db, Dp, cfm, roomwet1 + 1) + dxextra
                            ehr1 = eh - ((secload * 12000) / (latentfactoer * cfm))
                            sttemp1 = st(hre, elev, Dp, ehr1)
                            senstemp1 = db - ((secload * 12000) / (senfactoer * cfm))
                            if maxsen >= secload:
                                finltemp1 = senstemp1
                            else:
                                finltemp1 = sttemp1
                            finleh1 = ehr1
                            finlwet1 = weten(ehr1, elev)
                            finldry1 = finltemp1 - (finltemp1 - finlwet1) * deff

                            roomeh1 = finleh1 + ((roomload * 12000) / (latentfactoer * cfm))
                            roomwet1 = weten(roomeh1, elev)
                            roomdry1 = (
                                (roomload * 0.7 * 12000) / (senfactoer * cfm)
                            ) + finldry1

                            roomeh1 = round(roomeh1, 1)
                            roomwet1 = round(roomwet1, 1)
                            roomdry1 = round(roomdry1, 1)
                            if roomeh2 <= roomeh1:
                                break
                    else:
                        endtemp = roomdry * 3
                        return endtemp
                    if finldry1 > suptemp:
                        endtemp = roomdry1 * 5
                    else:
                        endtemp = roomdry1
                        return endtemp
                else:
                    endtemp = 1
                    return endtemp
            else:
                endtemp = 1
                return endtemp
        else:
            endtemp = 1
            return endtemp
    except :
        return DIV0


def rh(db, dew):
    try:
        db = float(db)
        dew = float(dew)
        es = 6.11 * 10 ** ((7.5 * ((db - 32) * 5 / 9)) / (237.7 + ((db - 32) * 5 / 9)))
        return (10 ** ((dew * 7.5) / (237.7 + dew))) * 611 / es
    except :
        return DIV0


def enthlpy(db, hr):
    db = float(db)
    hr = float(hr)
    fen = (0.24 * db) + (hr * (1061 + 0.444 * db))
    return fen


def enthlpywet(db, wb, elev):
    db = float(db)
    wb = float(wb)
    elev = float(elev)

    hr = humrat(db, wb, elev)
    return enthlpy(db, hr)


def weten(en, elev):
    try:
        en = float(en)
        elev = float(elev)

        dewa = 0.0003 * en ** 3 - 0.0491 * en ** 2 + 3.4088 * en - 1.6464 - 21

        for _ in range(1000):
            dewa = dewa + 0.01
            hra = (en - (0.24 * dewa)) / (1061 + 0.444 * dewa)
            dewa2 = dewpointal(elev, hra)
            if dewa >= dewa2:
                break
    except :
        return DIV0
    return dewa


def enddryroom(db, wb, elev, roomload, cfm, dxextra, deff):
    try:
        db = float(db)
        wb = float(wb)
        elev = float(elev)
        roomload = float(roomload)
        cfm = float(cfm)
        dxextra = float(dxextra)
        deff = float(deff)

        dbo = db
        senfactoer = 1.08 - (3.45454545454546e-05 * elev)
        latentfactoer = 4.45 - (0.00014234 * elev)
        rhe = relhum(db, wb, elev)
        Dp = dewpoint(db, rhe)
        hre = humrat(db, wb, elev)
        eh = enthlpywet(db, wb, elev)
        fstload = loadcalc(db, Dp, cfm, wb + 1) + dxextra
        ehr = eh - ((fstload * 12000) / (latentfactoer * cfm))
        maxsen = (senfactoer * cfm * (db - Dp)) / 12000
        sttemp = st(hre, elev, Dp, ehr)
        senstemp = db - ((fstload * 12000) / (senfactoer * cfm))

        if maxsen >= fstload:
            finltemp = senstemp
        else:
            finltemp = sttemp

        finlwet = weten(ehr, elev)
        finleh = ehr
        finldry = finltemp - (finltemp - finlwet) * deff
        roomeh = finleh + ((roomload * 12000) / (latentfactoer * cfm))
        roomwet = weten(roomeh, elev)
        roomdry = ((roomload * 0.7 * 12000) / (senfactoer * cfm)) + finldry
        roomeh1 = roomeh
        roomwet1 = roomwet
        roomdry1 = roomdry
        for _ in range(1000):
            roomeh2 = roomeh1
            secload = loadcalc(db, Dp, cfm, roomwet1 + 1) + dxextra
            ehr1 = eh - ((secload * 12000) / (latentfactoer * cfm))
            sttemp1 = st(hre, elev, Dp, ehr1)
            senstemp1 = db - ((secload * 12000) / (senfactoer * cfm))

            if maxsen >= secload:
                finltemp1 = senstemp1
            else:
                finltemp1 = sttemp1

            finleh1 = ehr1
            finlwet1 = weten(ehr1, elev)
            finldry1 = finltemp1 - (finltemp1 - finlwet1) * deff

            roomeh1 = finleh1 + ((roomload * 12000) / (latentfactoer * cfm))
            roomwet1 = weten(roomeh1, elev)
            roomdry1 = ((roomload * 0.7 * 12000) / (senfactoer * cfm)) + finldry1

            roomeh1 = round(roomeh1, 1)
            roomwet1 = round(roomwet1, 1)
            roomdry1 = round(roomdry1, 1)
            if roomeh2 <= roomeh1:
                break
        enddryroom = roomdry1
    except :
        return DIV0
    return enddryroom


def endwetroom(db, wb, elev, roomload, cfm, dxextra, deff):
    try:
        db = float(db)
        wb = float(wb)
        elev = float(elev)
        roomload = float(roomload)
        cfm = float(cfm)
        dxextra = float(dxextra)
        deff = float(deff)

        dbo = db
        senfactoer = 1.08 - (3.45454545454546e-05 * elev)
        latentfactoer = 4.45 - (0.00014234 * elev)
        rhe = relhum(db, wb, elev)
        Dp = dewpoint(db, rhe)
        hre = humrat(db, wb, elev)
        eh = enthlpywet(db, wb, elev)
        fstload = loadcalc(db, Dp, cfm, wb + 1) + dxextra
        ehr = eh - ((fstload * 12000) / (latentfactoer * cfm))
        maxsen = (senfactoer * cfm * (db - Dp)) / 12000
        sttemp = st(hre, elev, Dp, ehr)
        senstemp = db - ((fstload * 12000) / (senfactoer * cfm))

        if maxsen >= fstload:
            finltemp = senstemp
        else:
            finltemp = sttemp

        finlwet = weten(ehr, elev)
        finleh = ehr
        finldry = finltemp - (finltemp - finlwet) * deff
        roomeh = finleh + ((roomload * 12000) / (latentfactoer * cfm))
        roomwet = weten(roomeh, elev)
        roomdry = ((roomload * 0.7 * 12000) / (senfactoer * cfm)) + finldry
        roomeh1 = roomeh
        roomwet1 = roomwet
        roomdry1 = roomdry

        for _ in range(1000):
            roomeh2 = roomeh1
            roomwet2 = roomwet1
            roomdry2 = roomdry1
            secload = loadcalc(db, Dp, cfm, roomwet1 + 1) + dxextra
            ehr1 = eh - ((secload * 12000) / (latentfactoer * cfm))
            sttemp1 = st(hre, elev, Dp, ehr1)
            senstemp1 = db - ((secload * 12000) / (senfactoer * cfm))
            if maxsen >= secload:
                finltemp1 = senstemp1
            else:
                finltemp1 = sttemp1
            finleh1 = ehr1
            finlwet1 = weten(ehr1, elev)
            finldry1 = finltemp1 - (finltemp1 - finlwet1) * deff

            roomeh1 = finleh1 + ((roomload * 12000) / (latentfactoer * cfm))
            roomwet1 = weten(roomeh1, elev)
            roomdry1 = ((roomload * 0.7 * 12000) / (senfactoer * cfm)) + finldry1

            roomeh1 = round(roomeh1, 1)
            roomwet1 = round(roomwet1, 1)
            roomdry1 = round(roomdry1, 1)
            if roomeh2 <= roomeh1:
                break
        endwetroom = roomwet1
    except :
        return DIV0
    return endwetroom


def wetbulb(db, elev, rh):
    try:
        db = float(db)
        elev = float(elev)
        rh = float(rh)

        db = (db - 32) * 5 / 9
        asum1 = db
        wb = asum1
        pt = (14.696 * (1 - 0.0000068753 * elev) ** 5.2559) * 68.9476
        es = c1 * e ** ((c2 * db) / (db + c3))
        et1 = es * rh / 100
        et = et1
        for _ in range(1000):

            z1 = c1 * e ** ((c2 * wb) / (wb + c3))
            z2 = pt * (c4 + (c5 * wb)) * (db - wb)
            z3 = z1 - z2
            z = z3
            wb = wb - 1
            if z < et:
                break
        for _ in range(1000):
            y1 = c1 * e ** ((c2 * wb) / (wb + c3))
            y2 = pt * (c4 + (c5 * wb)) * (db - wb)
            y3 = y1 - y2
            y = y3
            wb = wb + 0.01
            if y >= et:
                break
        wb = (wb * 9 / 5) + 32
        wb = round(wb, 3)
        wetbulb = wb
    except :
        return DIV0
    return wetbulb


def findwin(db, Dp, cfm, ton):
    db = float(db)
    Dp = float(Dp)
    cfm = float(cfm)
    ton = float(ton)

    win = Dp - 30
    load1 = loadcalc(db, Dp, cfm, win)
    if load1 <= ton:
        return win
    for _ in range(1000):
        win = win + 0.1
        load1 = loadcalc(db, Dp, cfm, win)
        if load1 <= ton:
            break
    return win


def maxwet(dbr, dbmax, elev, eff):
    dbr = float(dbr)
    dbmax = float(dbmax)
    elev = float(elev)
    eff = float(eff)

    wmax = dbmax
    dewm = dbmax
    for _ in range(1000):
        rhmax = relhum(dbmax, wmax, elev)

        if rhmax >= 0.0264:
            dewm = dewpoint(dbmax, rhmax)
        else:
            dewm = 0
        dpi = dbmax - ((dbmax - wmax) * eff)
        rhi = rh(dpi, ((dewm - 32) * 5 / 9))
        Wiec = wetbulb(dpi, elev, rhi)
        if Wiec <= dbr:
            maxwet = wmax
            break
        else:
            wmax = wmax - 0.1
    return maxwet

def convert(x, u1: str, u2: str):
    x = float(x)
    if u1.lower() == 'm':
        return x * 3.28084
    elif u1.lower() == 'ft':
        return x / 3.28084
    elif u1.lower() == 'f':
        return (x - 32) * 5. / 9.
    elif u1.lower() == 'c':
        return (x * 9. / 5.) + 32
    else:
        return x