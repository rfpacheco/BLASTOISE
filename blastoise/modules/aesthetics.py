BOX_WIDTH_DEFAULT = 80


def print_message_box(message: str, width: int = BOX_WIDTH_DEFAULT) -> None:
    """
    Print ``message`` centered inside an ASCII art box.

    Parameters
    ----------
    message : str
        Text to display.
    width : int, optional
        Number of ``=`` characters in the top/bottom borders.
        If smaller than the message length, it is automatically enlarged.
    """
    width = max(width, len(message))  # ensure the message fits
    border = "=" * width

    print()  # leading blank line
    print(f"|{border}|")
    print(f"|{message:^{width}}|")
    print(f"|{border}|")
    print()  # trailing blank line


def blastoise_art() -> None:
    """
    Print ASCII art of Blastoise.
    """
    print("""
                       o O       o O       o O       o O       o O
                     o | | O   o | | O   o | | O   o | | O   o | | O
                   O | | | | O | | | | O | | | | O | | | | O | | | | O
                  O-oO | | o   O | | o   O | | o   O | | o   O | | oO-o
                 O---o O o       O o       O o       O o      O o  O---o
                O-----O                                           O-----o
                o-----O           ⣠⣴⣾⣶⣿⣿⣶⣶⣶⣿⡟⠀⠀⠀                  o-----O
                 o---O          ⣠⣼⣿⣿⡿⣋⣠⣿⣿⣿⣿⣿⡶⢶⣶⣤⣤⣀⣤⣶⣿⣗⡤⠶⢦⡀⠀        o---O
                  o-O          ⣤⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠀⠘⢿⣿⣿⣿⡿⢁⣾⡀⠀⢀⡷⠀         o-O
                   O           ⠻⣿⠟⠛⠛⠻⠟⠛⠋⢹⣿⣿⣿⢣⡆⠀⠈⠛⠛⠋⠀⢻⣿⣿⣶⠟⣻⣦          O  
                  o-O             ⣴⡿⠁   ⢸⣿⣿⢇⣿⠃⠀⣠⣤⣤⣤⣤⣀⠉⠙⢁⣴⡿⠁         o-O    
                 o---O       ⢀⣾⣿⣿ ⢠⠀⠀⠀⣰⣤⣿⣿⢋⣬⡄⢀⣾⣿⣿⣿⣿⣿⣿⣧⠀⣿⣯⠀⠀        o---O    
                O-----O  ⢀⣾⣿⣿⣿⣿   ⣠⠻⠿⠿⠿⠿⣛⣵⣿⣿⣧⢸⣿⣿⣿⣿⣿⣿⣿⣿⣄⣿⣿⡆⠀       O-----o    
                O-----O⢀⣾⣿⣿⣿⣿    ⣠⣿⡀⢸⣿⣿⣿⣿⣿⣿⣿⠿⠆⠻⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣷⠀       O-----o  
                O---⣿⣿⣿⣿⣿⣿⣿  ⢀⣀⣴⣾⣿⣿⡇⣬⣭⣭⣭⣭⣭⣶⣶⣿⣷⡄⢈⣻⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀       O-----o  
                o-⣿⣿⣿⣿⣿      ⠰⢾⣿⣿⣿⣿⡇⣿⣿⣿⣿⣿⣿⣿⣿⣿⡟⢐⣛⡻⣿⣿⣿⣿⣿⣿⠻⣿⣿⠀       o-----O
                 o-⣿⣿         ⠁⠀⣠⣶⣿⣷⢸⣿⣿⣿⣿⣿⣿⣿⡿⠿⠛⠋⡵⠿⢿⣿⣿⣿⢟⣄⢹⡏⠀        o---O
                  o-O          ⣰⣿⣿⣿⣿⣆⢲⣶⣶⣶⣶⣶⣶⣶⣿⢇⣷⣾⣿⡇⣟⣯⣶⣿⣿⡾⠀⠀         o-O  
                   O           ⣿⣿⣿⣿⣿⣿⣦⠹⣿⣿⣿⣿⣿⣿⣿⡜⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡾⠋          O
                  O-o          ⠘⣿⣿⣿⣿⣿⣿⣷⣬⠉⠿⣛⣻⣿⣯⣥⣹⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀         O-O
                 O---o        ⣠⣶⣿⣿⣿⣿⣿⣿⠿⠿⠦⠀⠀⠀⠉⠉⠁⠀⠹⣿⣿⣿⣿⣿⣿⣿⡿⠀⠀        O---o
                O-----o                         ⠀⠹⠿⠛⠿⣿⠟⠛⠛⠀        O-----o
                o-----O                                           o-----O
                 o---O o O       o O       o O       o O       o O o---O
                  o-Oo | | O   o | | O   o | | O   o | | O   o | | Oo-O
                   O | | | | O | | | | O | | | | O | | | | O | | | | O
                     O | | o   O | | o   O | | o   O | | o   O | | o
                       O o       O o       O o       O o       O o     
""")
