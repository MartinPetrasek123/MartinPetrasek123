document.getElementById('entryForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const location = document.getElementById('location').value;
    const coordinates = document.getElementById('coordinates').value;
    const note = document.getElementById('note').value;
    const photo = document.getElementById('photo').files[0];

    const reader = new FileReader();
    reader.onload = function () {
        const entry = {
            location,
            coordinates,
            note,
            photo: reader.result
        };
        const entries = JSON.parse(localStorage.getItem('entries') || '[]');
        entries.push(entry);
        localStorage.setItem('entries', JSON.stringify(entries));
        renderEntries();
    };

    if (photo) reader.readAsDataURL(photo);
    else reader.onload();
});

function renderEntries() {
    const container = document.getElementById('entries');
    container.innerHTML = '';
    const entries = JSON.parse(localStorage.getItem('entries') || '[]');
    entries.forEach(entry => {
        const div = document.createElement('div');
        div.className = 'entry';
        div.innerHTML = '<strong>' + entry.location + '</strong><br>' +
                        'Souřadnice: ' + entry.coordinates + '<br>' +
                        'Poznámka: ' + entry.note + '<br>' +
                        (entry.photo ? '<img src="' + entry.photo + '" width="200">' : '');
        container.appendChild(div);
    });
}
window.onload = renderEntries;