import krovak05


krovak = krovak05.Transformation()


with open("vstup.txt", "r") as f_in, open("vystup.txt", "w") as f_out:

    for l in f_in.readlines():

        x, y, h = map(lambda x: float(x), l.split(","))

        lat, lon, alt = krovak.jtsk_etrs(y, x, h)

        f_out.write(f"{lon:.15f},{lat:.15f},{alt:.5f}\n")
