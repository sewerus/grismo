from distutils.core import setup
import py2exe

setup(name="Grismo",
      version="0.1",
      author="Seweryn Panek",
      author_email="seweryn607@gmail.com",
      url="https://poehub.pl/",
      license="GNU General Public License (GPL)",
      packages=['grismo'],
      package_data={"grismo": ["ui/*"]},
      scripts=["bin/grismo"],
      windows=[{"script": "bin/grismo"}],
      options={"py2exe": {"skip_archive": True, "includes": ["sip"]}})